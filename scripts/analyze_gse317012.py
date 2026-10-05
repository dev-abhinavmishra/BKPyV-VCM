#!/usr/bin/env python3
"""Per-sample GSE317012 analysis and microarray series inventory.

Loads each sample's filtered 10x matrix with a sparse scipy path (no whole-
matrix densification), computes per-cell QC with transparent exclusion
accounting, marker-proxy cell-type fractions (clearly labelled as proxy
estimates), and five curated pathway scores. Sample-level summaries are
written to data/processed/, and condition contrasts are computed across
samples (the 26 biopsies are the replicates; cells are never treated as
independent patients). Also downloads/parses GEO series matrices for the
GSE47199/GSE75693/GSE72925 microarrays and emits an honest inventory of
what processed data are actually available.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
import tarfile
import time
import tracemalloc
import urllib.request
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scipy.stats

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = REPO_ROOT / "data" / "raw" / "GSE317012_RAW"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
RESEARCH_DIR = REPO_ROOT / "data" / "research"
GSE317012_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317012/matrix/"
    "GSE317012_series_matrix.txt.gz"
)

MIN_GENES_PER_CELL = 200
MAX_MT_FRACTION = 0.20
NORMALIZATION_TARGET = 1e4

MISSING_RAW_NOTE = """\
ERROR: required raw input is missing: {path}

The raw 10x extraction tree (data/raw/GSE317012_RAW/) and the source archive
(data/GSE317012_RAW.tar, sha256
1EA67E31A190AA05F39A7664721BC0CAF3D219997DFE67C2FC7920C4854327A8) were
deleted on 2026-09-19 under user-authorized storage reclamation. The
processed outputs and provenance records under data/processed/ are the
retained evidence.

To restore and rerun:
  1. Re-download the archive from NCBI GEO FTP:
     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317012/suppl/GSE317012_RAW.tar
     (verify sha256 == 1EA67E31...4327A8, recorded in
     data/processed/gse317012_provenance.csv)
  2. Extract + organize:
     python scripts/reorganize_gse317012.py data/raw/GSE317012_RAW \\
       --extract data/GSE317012_RAW.tar --move --verify
  3. Re-run this analyzer.
"""

# Pathway gene sets curated in docs/bkpyv_research_to_model_map.md
# (S/G2M, DDR, mitochondrial) plus canonical sets for the pathways the map
# names but leaves unenumerated (translation, antigen presentation).
GENE_SETS: dict[str, list] = {
    "s_g2m": ["CLSPN", "TOP2A", "MKI67"],
    "ddr": ["BRCA1", "BRCA2", "PRKDC", "FANCI", "MMS22L"],
    "mito": ["MT-ND4", "MT-CO1", "MT-CYB", "MT-ATP6"],
    "translation": [
        "RPL3",
        "RPL5",
        "RPL7",
        "RPL8",
        "RPLP0",
        "RPLP1",
        "RPS3",
        "RPS6",
        "RPS18",
        "RPS27A",
    ],
    "antigen_presentation": [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "B2M",
        "TAP1",
        "TAP2",
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DPA1",
    ],
}
MARKER_SETS: dict[str, list] = {
    "epithelial_proxy": ["LRP2", "CUBN", "UMOD", "SLC34A1"],
    "immune_proxy": ["PTPRC"],
}

MICROARRAY_SERIES = {
    "GSE47199": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE471nnn/GSE47199/matrix/GSE47199_series_matrix.txt.gz",
    "GSE75693": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75693/matrix/GSE75693_series_matrix.txt.gz",
    "GSE72925": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72925/matrix/GSE72925_series_matrix.txt.gz",
}

TISSUE_KEYWORDS = ("blood", "biopsy", "urine", "kidney", "renal")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _open_maybe_gz(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _unquote(field: str) -> str:
    field = field.strip()
    if len(field) >= 2 and field.startswith('"') and field.endswith('"'):
        return field[1:-1]
    return field


def download(url: str, dest: Path) -> dict:
    """Download url to dest (no overwrite of existing different content)."""
    fetched = datetime.now(timezone.utc).isoformat(timespec="seconds")
    req = urllib.request.Request(url, headers={"User-Agent": "vcm-analysis/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = resp.read()
    if dest.exists() and dest.read_bytes() != payload:
        raise ValueError(f"refusing to overwrite existing file with different content: {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        dest.write_bytes(payload)
    return {
        "url": url,
        "fetched_utc": fetched,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "bytes": len(payload),
        "path": str(dest.relative_to(REPO_ROOT)),
    }


# ---------- series matrix parsing ----------


def parse_series_matrix(path: Path) -> dict[str, dict]:
    """Parse a GEO series-matrix file into {GSM: sample metadata dict}."""
    fields: dict[str, list] = {}
    has_values = False
    with _open_maybe_gz(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!Series_matrix_table_begin"):
                has_values = True
                continue
            if not line.startswith("!Sample_"):
                continue
            parts = [_unquote(p) for p in line.split("\t")]
            fields.setdefault(parts[0], []).append(parts[1:])

    accessions = fields.get("!Sample_geo_accession", [[]])[0]
    samples: dict[str, dict] = {}
    for i, gsm in enumerate(accessions):
        characteristics = [
            row[i] for row in fields.get("!Sample_characteristics_ch1", []) if i < len(row)
        ]
        samples[gsm] = {
            "gsm": gsm,
            "title": (
                (fields.get("!Sample_title", [[]])[0] + [""] * (i + 1))[i]
                if fields.get("!Sample_title")
                else ""
            ),
            "source_name": (
                (fields.get("!Sample_source_name_ch1", [[]])[0] + [""] * (i + 1))[i]
                if fields.get("!Sample_source_name_ch1")
                else ""
            ),
            "characteristics": characteristics,
            "phase": derive_phase(characteristics),
        }
    for s in samples.values():
        s["has_expression_values"] = has_values
    return samples


def derive_phase(characteristics: Iterable[str]) -> str:
    """Derive the disease phase from authoritative characteristics text."""
    text = " | ".join(characteristics).lower()
    for phase in ("control", "peaking", "resolving"):
        if phase in text:
            return phase.capitalize()
    return "Unclassified"


def crosscheck_gsm_sets(matrix_gsms: set, tar_gsms: set) -> None:
    """Require the metadata GSM set and the raw-data GSM set to agree."""
    if matrix_gsms == tar_gsms:
        return
    only_matrix = sorted(matrix_gsms - tar_gsms)
    only_tar = sorted(tar_gsms - matrix_gsms)
    raise ValueError(
        "GSM set mismatch between series matrix and raw archive: "
        f"matrix-only={only_matrix} tar-only={only_tar}"
    )


# ---------- sparse 10x loading / QC / scoring ----------


def load_filtered_10x(sample_dir: Path) -> tuple:
    """Load (matrix genes x_counts cells csr, feature names, barcodes)."""
    matrix = None
    for name in ("matrix.mtx.gz", "matrix.mtx"):
        p = sample_dir / name
        if p.is_file():
            src = gzip.open(p, "rb") if name.endswith(".gz") else open(p, "rb")
            with src:
                matrix = scipy.io.mmread(src).tocsr()
            break
    if matrix is None:
        raise FileNotFoundError(f"no filtered matrix.mtx in {sample_dir}")

    def read_lines(prefix):
        for name in (f"{prefix}.tsv.gz", f"{prefix}.tsv"):
            p = sample_dir / name
            if p.is_file():
                with _open_maybe_gz(p) as fh:
                    return [ln.rstrip("\n") for ln in fh]
        raise FileNotFoundError(f"no {prefix}.tsv in {sample_dir}")

    feature_rows = read_lines("features")
    features = [r.split("\t")[1] if "\t" in r else r.split("\t")[0] for r in feature_rows]
    barcodes = read_lines("barcodes")

    if matrix.shape[0] == len(barcodes) and matrix.shape[1] == len(features):
        matrix = matrix.T.tocsr()
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            f"matrix shape {matrix.shape} incompatible with {len(features)} features / "
            f"{len(barcodes)} barcodes in {sample_dir}"
        )
    return matrix, features, barcodes


def qc_filter(
    x_counts: sp.csr_matrix,
    features: list,
    min_genes: int = MIN_GENES_PER_CELL,
    max_mt_frac: float = MAX_MT_FRACTION,
) -> tuple:
    """Per-cell QC on genes x_counts cells csr. Returns (X_keep, stats)."""
    n_genes = np.asarray((x_counts > 0).sum(axis=0)).ravel()
    counts = np.asarray(x_counts.sum(axis=0)).ravel()
    mt_idx = [i for i, g in enumerate(features) if g.startswith("MT-")]
    mt_counts = (
        np.asarray(x_counts[mt_idx, :].sum(axis=0)).ravel()
        if mt_idx
        else np.zeros(x_counts.shape[1])
    )
    mt_frac = np.divide(
        mt_counts, counts, out=np.zeros(x_counts.shape[1], dtype=float), where=counts > 0
    )
    keep = (n_genes >= min_genes) & (mt_frac <= max_mt_frac)
    stats = {
        "cells_pre": x_counts.shape[1],
        "cells_post": int(keep.sum()),
        "excluded": int((~keep).sum()),
        "excluded_low_genes": int((n_genes < min_genes).sum()),
        "excluded_high_mt": int((mt_frac > max_mt_frac).sum()),
    }
    return x_counts[:, keep].tocsr(), stats


def normalize_log1p(x_counts: sp.csr_matrix) -> sp.csr_matrix:
    """Library-size normalize each cell to 1e4 then log1p (stays sparse)."""
    lib = np.asarray(x_counts.sum(axis=0)).ravel()
    scale = np.divide(NORMALIZATION_TARGET, lib, out=np.zeros(lib.size, dtype=float), where=lib > 0)
    x_norm = x_counts.multiply(sp.csr_matrix(scale)).tocsr()
    x_norm.data = np.log1p(x_norm.data)
    return x_norm


def gene_set_scores(
    x_norm: sp.csr_matrix, features: list, gene_sets: Mapping[str, list]
) -> dict[str, np.ndarray]:
    """Per-cell score = mean normalized expression over the set's present genes."""
    name_to_row = {g: i for i, g in enumerate(features)}
    scores = {}
    for name, genes in gene_sets.items():
        rows = [name_to_row[g] for g in genes if g in name_to_row]
        if not rows:
            scores[name] = np.full(x_norm.shape[1], np.nan)
            continue
        scores[name] = np.asarray(x_norm[rows, :].mean(axis=0)).ravel()
    return scores


def marker_fractions(
    x_counts: sp.csr_matrix, features: list, marker_sets: Mapping[str, list]
) -> dict[str, float]:
    """Fraction of cells expressing >=1 marker gene (raw counts > 0)."""
    name_to_row = {g: i for i, g in enumerate(features)}
    n_cells = x_counts.shape[1]
    fracs = {}
    for name, genes in marker_sets.items():
        rows = [name_to_row[g] for g in genes if g in name_to_row]
        if not rows or n_cells == 0:
            fracs[name] = np.nan
            continue
        expressing = np.asarray((x_counts[rows, :] > 0).sum(axis=0)).ravel() > 0
        fracs[name] = float(expressing.sum() / n_cells)
    return fracs


def analyze_sample(sample_dir: Path, meta: Mapping) -> dict:
    """One GSM: load filtered 10x, QC, fractions, pathway scores."""
    x_counts, features, barcodes = load_filtered_10x(sample_dir)
    x_qc, stats = qc_filter(x_counts, features)
    x_norm = normalize_log1p(x_qc)
    scores = gene_set_scores(x_norm, features, GENE_SETS)
    fracs = marker_fractions(x_qc, features, MARKER_SETS)
    curated = set(g for gs in list(GENE_SETS.values()) + list(MARKER_SETS.values()) for g in gs)
    missing = sorted(curated - set(features))
    row = {
        "gsm": meta.get("gsm", sample_dir.name),
        "phase": meta.get("phase", "Unclassified"),
        "title": meta.get("title", ""),
        "source_name": meta.get("source_name", ""),
        "cells_pre_qc": stats["cells_pre"],
        "cells_post_qc": stats["cells_post"],
        "cells_excluded": stats["excluded"],
        "barcodes_in_file": len(barcodes),
        "frac_epithelial_proxy": round(fracs["epithelial_proxy"], 6),
        "frac_immune_proxy": round(fracs["immune_proxy"], 6),
        "missing_markers": ";".join(missing),
    }
    for name, per_cell in scores.items():
        valid = per_cell[~np.isnan(per_cell)]
        row[f"score_{name}"] = round(float(valid.mean()), 6) if valid.size else np.nan
    return row


# ---------- contrasts ----------


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """P(a>b) - P(a<b); distribution-free effect size."""
    gt = sum(x_counts > y for x_counts in a for y in b)
    lt = sum(x_counts < y for x_counts in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def contrast(df: pd.DataFrame, column: str, group_a: str, group_b: str) -> dict:
    a = df.loc[df["phase"] == group_a, column].to_numpy(dtype=float)
    b = df.loc[df["phase"] == group_b, column].to_numpy(dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    u_stat, p = scipy.stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "variable": column,
        "group_a": group_a,
        "group_b": group_b,
        "n_a": len(a),
        "n_b": len(b),
        "mean_a": round(float(a.mean()), 6),
        "mean_b": round(float(b.mean()), 6),
        "mean_diff_a_minus_b": round(float(a.mean() - b.mean()), 6),
        "cliffs_delta": round(cliffs_delta(a, b), 6),
        "mwu_stat": round(float(u_stat), 6),
        "mwu_p": p,
    }


def run_contrasts(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in df.columns if c.startswith(("score_", "frac_"))]
    rows = []
    for col in cols:
        for a, b in (("Peaking", "Control"), ("Resolving", "Peaking")):
            rows.append(contrast(df, col, a, b))
    return pd.DataFrame(rows)


# ---------- microarray inventory ----------


def series_matrix_url(gse: str) -> str:
    prefix = f"{gse[:-3]}nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/"
        f"{gse}_series_matrix.txt.gz"
    )


def tissue_context(samples: Mapping[str, dict]) -> str:
    seen = set()
    for s in samples.values():
        text = f"{s.get('source_name', '')} {' '.join(s.get('characteristics', []))}".lower()
        for kw in TISSUE_KEYWORDS:
            if kw in text:
                seen.add(kw)
    return ";".join(sorted(seen)) if seen else "unknown"


def tar_gsm_members(tar_path: Path) -> list:
    if not tar_path.is_file():
        return []
    with tarfile.open(tar_path) as tf:
        return [m.name for m in tf.getmembers() if m.isfile()]


def inventory_series(gse: str) -> dict:
    """Fetch series matrix, summarize contents, and check local raw tar."""
    url = series_matrix_url(gse)
    dest = RESEARCH_DIR / f"{gse}_series_matrix.txt.gz"
    tar = REPO_ROOT / "data" / f"{gse}_RAW.tar"
    row = {
        "series": gse,
        "matrix_url": url,
        "matrix_fetched_utc": "",
        "matrix_sha256": "",
        "matrix_status": "unreachable",
        "n_gsm_in_matrix": 0,
        "n_files_in_local_tar": len(tar_gsm_members(tar)),
        "platform": "",
        "tissue_context": "",
        "has_expression_values": False,
    }
    try:
        meta = download(url, dest)
    except urllib.error.HTTPError as exc:
        row["matrix_status"] = f"HTTP {exc.code}"
        return row
    except (urllib.error.URLError, OSError, TimeoutError) as exc:
        row["matrix_status"] = f"fetch-error {type(exc).__name__}: {exc}"
        return row
    row.update(
        {
            "matrix_fetched_utc": meta["fetched_utc"],
            "matrix_sha256": meta["sha256"],
            "matrix_status": "downloaded",
        }
    )
    samples = parse_series_matrix(dest)
    row["n_gsm_in_matrix"] = len(samples)
    row["tissue_context"] = tissue_context(samples)
    row["has_expression_values"] = any(s["has_expression_values"] for s in samples.values())
    platform_fields = []
    with _open_maybe_gz(dest) as fh:
        for line in fh:
            if line.startswith("!Series_platform_id"):
                platform_fields = [_unquote(p) for p in line.rstrip("\n").split("\t")[1:]]
                break
    row["platform"] = platform_fields[0] if platform_fields else ""
    summary = pd.DataFrame(
        [
            {
                "gsm": s["gsm"],
                "title": s["title"],
                "source_name": s["source_name"],
                "characteristics": " | ".join(s["characteristics"]),
                "platform": row["platform"],
            }
            for s in samples.values()
        ]
    )
    summary.to_csv(PROCESSED_DIR / f"{gse.lower()}_series_summary.csv", index=False)
    return row


# ---------- main pipeline ----------


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, default=RAW_ROOT)
    parser.add_argument("--processed-dir", type=Path, default=PROCESSED_DIR)
    parser.add_argument("--skip-microarray-inventory", action="store_true")
    args = parser.parse_args(argv)
    args.processed_dir.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    tracemalloc.start()
    downloads = {}

    meta_dest = RESEARCH_DIR / "gse317012_series_matrix.txt.gz"
    meta = download(GSE317012_MATRIX_URL, meta_dest)
    downloads[meta["path"]] = meta
    samples_meta = parse_series_matrix(meta_dest)

    if not args.raw_root.is_dir() or not any(
        p.is_dir() and p.name.startswith("GSM") for p in args.raw_root.iterdir()
    ):
        print(MISSING_RAW_NOTE.format(path=args.raw_root), file=sys.stderr)
        return 1

    sample_dirs = sorted(
        p for p in args.raw_root.iterdir() if p.is_dir() and p.name.startswith("GSM")
    )
    tar_gsms = {p.name for p in sample_dirs}
    crosscheck_gsm_sets(set(samples_meta), tar_gsms)

    meta_rows = [
        {
            "gsm": s["gsm"],
            "title": s["title"],
            "phase": s["phase"],
            "source_name": s["source_name"],
            "characteristics": " | ".join(s["characteristics"]),
        }
        for s in (samples_meta[g] for g in sorted(tar_gsms))
    ]
    meta_df = pd.DataFrame(meta_rows)
    meta_csv = args.processed_dir / "gse317012_sample_metadata.csv"
    meta_df.to_csv(meta_csv, index=False)
    phase_counts = meta_df["phase"].value_counts().to_dict()

    rows, total_cells = [], 0
    for sd in sample_dirs:
        row = analyze_sample(sd, samples_meta[sd.name])
        total_cells += row["cells_post_qc"]
        rows.append(row)
    scores_df = pd.DataFrame(rows)
    scores_csv = args.processed_dir / "gse317012_pathway_scores.csv"
    scores_df.to_csv(scores_csv, index=False)

    contrasts_df = run_contrasts(scores_df)
    contrasts_csv = args.processed_dir / "gse317012_contrasts.csv"
    contrasts_df.to_csv(contrasts_csv, index=False)

    provenance_rows = [
        {
            "file": meta["path"],
            "sha256": meta["sha256"],
            "role": "input-download",
            "detail": meta["url"] + " @ " + meta["fetched_utc"],
        },
    ]
    for sd in sample_dirs:
        for f in sorted(sd.iterdir()):
            provenance_rows.append(
                {
                    "file": str(f.relative_to(REPO_ROOT)),
                    "sha256": sha256_file(f),
                    "role": "input-10x",
                    "detail": sd.name,
                }
            )
    for csv_path in (meta_csv, scores_csv, contrasts_csv):
        provenance_rows.append(
            {
                "file": str(csv_path.relative_to(REPO_ROOT)),
                "sha256": sha256_file(csv_path),
                "role": "output",
                "detail": "generated by scripts/analyze_gse317012.py",
            }
        )

    inventory_rows = []
    if not args.skip_microarray_inventory:
        for gse in MICROARRAY_SERIES:
            inv = inventory_series(gse)
            inventory_rows.append(inv)
            if inv["matrix_sha256"]:
                dl = RESEARCH_DIR / f"{gse}_series_matrix.txt.gz"
                provenance_rows.append(
                    {
                        "file": str(dl.relative_to(REPO_ROOT)),
                        "sha256": inv["matrix_sha256"],
                        "role": "input-download",
                        "detail": inv["matrix_url"] + " @ " + inv["matrix_fetched_utc"],
                    }
                )
        inv_csv = args.processed_dir / "microarray_inventory.csv"
        pd.DataFrame(inventory_rows).to_csv(inv_csv, index=False)
        provenance_rows.append(
            {
                "file": str(inv_csv.relative_to(REPO_ROOT)),
                "sha256": sha256_file(inv_csv),
                "role": "output",
                "detail": "microarray series inventory",
            }
        )
        for gse in MICROARRAY_SERIES:
            summary = args.processed_dir / f"{gse.lower()}_series_summary.csv"
            if summary.is_file():
                provenance_rows.append(
                    {
                        "file": str(summary.relative_to(REPO_ROOT)),
                        "sha256": sha256_file(summary),
                        "role": "output",
                        "detail": f"{gse} per-GSM series summary",
                    }
                )

    prov_csv = args.processed_dir / "gse317012_provenance.csv"
    pd.DataFrame(provenance_rows).to_csv(prov_csv, index=False)

    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    runtime = time.perf_counter() - started
    print("phase_counts:", phase_counts)
    print(f"samples_analyzed: {len(rows)}")
    print(f"cells_analyzed_post_qc: {total_cells}")
    print(f"runtime_s: {runtime:.1f}")
    print(f"peak_python_alloc_mb: {peak / 1e6:.1f}")
    print(f"wrote {scores_csv}")
    print(f"wrote {contrasts_csv}")
    print(f"wrote {prov_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
