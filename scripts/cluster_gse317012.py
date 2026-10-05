#!/usr/bin/env python3
"""Cell-level clustering, annotation, and infected-vs-bystander analysis for GSE317012.

Complements ``scripts/analyze_gse317012.py`` (which works at biopsy level) with a
single-cell pipeline over the same QC-filtered 10x data:

  QC -> log-normalization -> HVG -> PCA -> neighbors -> UMAP -> Leiden ->
  marker annotation -> viral-response module score -> infected-signature
  vs bystander comparison within tubular epithelial cells.

IMPORTANT INTERPRETATION BOUNDARY
---------------------------------
GSE317012 was quantified against a human-only reference (36,601 Ensembl genes;
no BKPyV genes are present in ``features.tsv.gz``, verified 2026-10-05). Viral
transcripts therefore cannot be counted directly. "Infected-signature" cells
are tubular epithelial cells in peaking/resolving biopsies whose composite
host viral-response module score falls in the top decile of the
phase-matched epithelial distribution — a *signature proxy*, not a direct
viral measurement. Control-biopsy epithelial cells are the reference
population for "uninfected".

Outputs (all under data/processed/ unless noted):

  gse317012_cells.parquet            per-cell: gsm, phase, leiden, cell_type,
                                     module scores, infected_call
  gse317012_cluster_summary.csv      per-cluster counts + marker means + label
  gse317012_infected_vs_bystander_de.csv  Wilcoxon DE, infected vs bystander
  gse317012_model_mapping.csv        finding -> ODE state variable/parameter
  gse317012_calibration_holdout.csv  which GSMs are calibration vs validation
  gse317012_calibration.csv          pathway-score ratios -> suggested ODE params
  outputs/figures/single_cell/*.png  UMAP / score / heatmap figures

Deterministic: fixed seed (42), deterministic holdout hash-split.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = REPO_ROOT / "data" / "raw" / "GSE317012_RAW"
RESEARCH_DIR = REPO_ROOT / "data" / "research"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
FIGURE_DIR = REPO_ROOT / "outputs" / "figures" / "single_cell"
SERIES_MATRIX = RESEARCH_DIR / "gse317012_series_matrix.txt.gz"

sys.path.insert(0, str(REPO_ROOT / "scripts"))
from analyze_gse317012 import (  # noqa: E402
    GENE_SETS,
    MIN_GENES_PER_CELL,
    MAX_MT_FRACTION,
    NORMALIZATION_TARGET,
    parse_series_matrix,
)

SEED = 42

# Gene sets added for the infected-vs-bystander readout (paper summary: immune
# signaling, wound healing, ECM remodeling, metabolic restructuring).
EXTRA_GENE_SETS: dict[str, list[str]] = {
    "interferon": [
        "MX1",
        "ISG15",
        "OAS1",
        "OAS2",
        "OAS3",
        "IFI6",
        "IFIT1",
        "IFIT3",
        "IRF7",
        "STAT1",
    ],
    "ecm_wound": [
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "FN1",
        "VIM",
        "TGFB1",
        "TIMP1",
        "MMP7",
    ],
}

# Canonical kidney cell-type marker panels for cluster annotation.
CELL_TYPE_MARKERS: dict[str, list[str]] = {
    "proximal_tubule": ["LRP2", "CUBN", "SLC34A1", "SLC5A2", "SLC22A6"],
    "loop_of_henle": ["UMOD", "SLC12A1", "CLDN16"],
    "distal_convoluted_tubule": ["SLC12A3", "SLC8A1"],
    "collecting_duct": ["AQP2", "AQP3", "SLC4A1"],
    "podocyte": ["NPHS1", "NPHS2", "PODXL"],
    "endothelial": ["PECAM1", "KDR", "VWF"],
    "fibroblast": ["COL1A1", "DCN", "LUM"],
    "t_cell": ["CD3D", "CD3E", "CD8A", "CD4"],
    "nk_cell": ["NKG7", "GNLY", "KLRD1"],
    "b_cell": ["CD79A", "MS4A1", "CD19"],
    "myeloid": ["LST1", "FCGR3A", "S100A8", "LYZ"],
    "immune_other": ["PTPRC"],
}

TUBULAR_TYPES = {
    "proximal_tubule",
    "loop_of_henle",
    "distal_convoluted_tubule",
    "collecting_duct",
}


def load_and_qc() -> "tuple":
    """Load all GSM filtered 10x matrices, apply QC identical to the
    biopsy-level analyzer (>=200 genes, <=20% mito), tag cells with
    gsm + phase, and return the concatenated AnnData."""
    import scanpy as sc

    samples_meta = parse_series_matrix(SERIES_MATRIX)
    sample_dirs = sorted(
        p for p in RAW_ROOT.iterdir() if p.is_dir() and p.name.startswith("GSM")
    )

    adatas = []
    for sd in sample_dirs:
        adata = sc.read_10x_mtx(str(sd), var_names="gene_symbols")
        adata.var_names_make_unique()
        adata.obs["gsm"] = sd.name
        adata.obs["phase"] = samples_meta[sd.name]["phase"]
        adatas.append(adata)

    import anndata

    adata = anndata.concat(
        adatas, join="outer", index_unique=None, fill_value=0
    )
    del adatas

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    n_pre = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= MIN_GENES_PER_CELL)
        & (adata.obs["pct_counts_mt"] <= MAX_MT_FRACTION * 100.0)
    ].copy()
    print(f"QC: {n_pre} -> {adata.n_obs} cells "
          f"(min_genes={MIN_GENES_PER_CELL}, max_mt={MAX_MT_FRACTION})")
    return adata


def cluster_and_annotate(adata):
    """Normalize, embed, cluster, and annotate cell types by marker panels."""
    import scanpy as sc

    sc.pp.normalize_total(adata, target_sum=NORMALIZATION_TARGET)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep normalized values for scoring/DE

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    # PCA on the HVG subset only (sparse, non-centred SVD) — keeps every gene
    # available downstream for module scoring and DE while avoiding a dense
    # 35k x 36.6k scaled matrix.
    sc.pp.pca(
        adata,
        n_comps=50,
        use_highly_variable=True,
        zero_center=False,
        random_state=SEED,
    )
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED, flavor="igraph",
                 n_iterations=2, directed=False)

    # Marker annotation: score each cluster by mean normalized expression
    # of each cell-type panel; assign the best-scoring type.
    cluster_scores = {}
    for ctype, genes in CELL_TYPE_MARKERS.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        sc.tl.score_genes(
            adata, present, score_name=f"mk_{ctype}", random_state=SEED
        )
        cluster_scores[ctype] = f"mk_{ctype}"

    per_cluster = (
        adata.obs.groupby("leiden")[[c for c in cluster_scores.values()]]
        .mean()
    )
    labels = {}
    for cl, row in per_cluster.iterrows():
        best = row.idxmax()
        labels[cl] = best.removeprefix("mk_")
    adata.obs["cell_type"] = adata.obs["leiden"].map(labels)
    return adata, per_cluster


def score_modules(adata):
    """Score every cell for the curated pathway modules."""
    import scanpy as sc

    all_sets = {**GENE_SETS, **EXTRA_GENE_SETS}
    for name, genes in all_sets.items():
        present = [g for g in genes if g in adata.var_names]
        missing = sorted(set(genes) - set(present))
        if missing:
            print(f"  module {name}: missing genes {missing}")
        if not present:
            continue
        sc.tl.score_genes(
            adata, present, score_name=f"mod_{name}", random_state=SEED
        )
    # Composite viral-response score: hallmarks named by the source paper —
    # translation/ribosome biogenesis, cell cycle, DDR, interferon,
    # antigen presentation. Mito is excluded from the composite (it tracks
    # viability as much as infection).
    composite_parts = [
        c
        for c in (
            "mod_translation",
            "mod_s_g2m",
            "mod_ddr",
            "mod_interferon",
            "mod_antigen_presentation",
        )
        if c in adata.obs.columns
    ]
    adata.obs["viral_response_score"] = adata.obs[composite_parts].mean(axis=1)
    return adata


def call_infected(adata) -> pd.DataFrame:
    """Label infected-signature vs bystander tubular epithelial cells.

    Only tubular-epithelial clusters in peaking/resolving biopsies are
    eligible; the top decile of the phase-matched epithelial viral-response
    score is called 'infected_signature', the rest 'bystander'. Control
    epithelial cells are a separate reference class.
    """
    obs = adata.obs
    epi = obs["cell_type"].isin(TUBULAR_TYPES)
    infected_phase = obs["phase"].isin(["Peaking", "Resolving"])
    eligible = epi & infected_phase

    cutoff = obs.loc[eligible, "viral_response_score"].quantile(0.90)
    call = pd.Series("other", index=obs.index, dtype=object)
    call[epi & ~infected_phase] = "control_epithelial"
    call[eligible & (obs["viral_response_score"] >= cutoff)] = (
        "infected_signature"
    )
    call[eligible & (obs["viral_response_score"] < cutoff)] = "bystander"
    obs["infected_call"] = call
    print(
        "infected_call counts:",
        call.value_counts().to_dict(),
        f"(cutoff={cutoff:.3f})",
    )
    return obs


def run_de(adata) -> pd.DataFrame:
    """Wilcoxon DE: infected-signature vs bystander tubular epithelial."""
    import scanpy as sc

    mask = adata.obs["infected_call"].isin(["infected_signature", "bystander"])
    sub = adata[mask].copy()
    sc.tl.rank_genes_groups(
        sub,
        "infected_call",
        groups=["infected_signature"],
        reference="bystander",
        method="wilcoxon",
        use_raw=True,
    )
    de = sc.get.rank_genes_groups_df(sub, group="infected_signature")
    return de


def split_samples(samples_meta: dict) -> pd.DataFrame:
    """Deterministic stratified holdout: ~1/3 of GSMs per phase are
    validation, the rest calibration. Sorted + seeded so it reproduces."""
    rng = np.random.default_rng(SEED)
    rows = []
    for phase in sorted({s["phase"] for s in samples_meta.values()}):
        gsms = sorted(g for g, s in samples_meta.items() if s["phase"] == phase)
        order = rng.permutation(len(gsms))
        n_val = max(1, round(len(gsms) / 3))
        val_idx = set(order[:n_val])
        for i, gsm in enumerate(gsms):
            rows.append(
                {
                    "gsm": gsm,
                    "phase": phase,
                    "split": "validation" if i in val_idx else "calibration",
                }
            )
    return pd.DataFrame(rows)


def calibrate(obs: pd.DataFrame, splits: pd.DataFrame) -> pd.DataFrame:
    """Map measured module-score shifts onto ODE parameters/state variables.

    Only calibration-split samples feed the suggested values; the validation
    split is reported separately so the reader can check direction agreement.
    """
    modules = {
        "mod_translation": "translation_enhancement (ODE param; pathway_translation state)",
        "mod_ddr": "ddr_enhancement (ODE param; pathway_dna_damage_response state)",
        "mod_mito": "mitochondrial_importance (ODE param; pathway_mitochondrial_stress state)",
        "mod_interferon": "p_immune_prod / IFN (ODE params/state; pathway_interferon_response)",
        "mod_antigen_presentation": "pathway_antigen_presentation state",
        "mod_ecm_wound": "pathway_cellular_stress state (wound/ECM proxy)",
        "mod_s_g2m": "s_phase_bonus / CC gate (ODE params; pathway_cell_cycle)",
    }
    calib_gsms = set(splits.loc[splits["split"] == "calibration", "gsm"])
    valid_gsms = set(splits.loc[splits["split"] == "validation", "gsm"])

    rows = []
    epi = obs["infected_call"].isin(
        ["infected_signature", "bystander", "control_epithelial"]
    )
    for mod, target in modules.items():
        for split_name, gsm_set in (
            ("calibration", calib_gsms),
            ("validation", valid_gsms),
        ):
            sub = obs[epi & obs["gsm"].isin(gsm_set)]
            inf = sub.loc[sub["infected_call"] == "infected_signature", mod]
            ctl = sub.loc[sub["infected_call"] == "control_epithelial", mod]
            if len(inf) < 50 or len(ctl) < 50:
                continue
            rows.append(
                {
                    "module": mod,
                    "model_target": target,
                    "split": split_name,
                    "mean_infected_signature": float(inf.mean()),
                    "mean_control_epithelial": float(ctl.mean()),
                    "shift": float(inf.mean() - ctl.mean()),
                    "n_infected": int(len(inf)),
                    "n_control": int(len(ctl)),
                }
            )
    df = pd.DataFrame(rows)
    # Suggested parameter direction: qualitative only — a score shift is not
    # a rate constant; the mapping table records the evidence.
    return df


def make_figures(adata, out_dir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import scanpy as sc

    sc.settings.figdir = out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    sc.pl.umap(adata, color=["cell_type"], save="_celltype.png", show=False)
    sc.pl.umap(adata, color=["phase"], save="_phase.png", show=False)
    sc.pl.umap(
        adata, color=["viral_response_score"], save="_viral_score.png",
        show=False, cmap="viridis",
    )
    sc.pl.violin(
        adata[adata.obs["cell_type"].isin(TUBULAR_TYPES)],
        keys=["viral_response_score"],
        groupby="phase",
        save="_epithelial_score_by_phase.png",
        show=False,
        rotation=45,
    )


def main(argv: list[str] | None = None) -> int:
    global RAW_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, default=RAW_ROOT)
    parser.add_argument("--processed-dir", type=Path, default=PROCESSED_DIR)
    parser.add_argument("--figure-dir", type=Path, default=FIGURE_DIR)
    args = parser.parse_args(argv)

    if not args.raw_root.is_dir():
        print(
            f"ERROR: raw 10x tree missing: {args.raw_root}. "
            "Run scripts/download_gse317012.py first.",
            file=sys.stderr,
        )
        return 1
    RAW_ROOT = args.raw_root

    args.processed_dir.mkdir(parents=True, exist_ok=True)

    adata = load_and_qc()
    adata, per_cluster = cluster_and_annotate(adata)
    adata = score_modules(adata)
    obs = run_holdout_and_calls(adata, args.processed_dir)
    de = run_de(adata)
    de.to_csv(args.processed_dir / "gse317012_infected_vs_bystander_de.csv",
              index=False)
    make_figures(adata, args.figure_dir)
    print("done")
    return 0


def run_holdout_and_calls(adata, processed_dir: Path):
    """Call infected vs bystander, split samples, calibrate, write tables."""
    obs = call_infected(adata)

    samples_meta = parse_series_matrix(SERIES_MATRIX)
    splits = split_samples(samples_meta)
    splits.to_csv(processed_dir / "gse317012_calibration_holdout.csv",
                  index=False)

    calib = calibrate(obs, splits)
    calib.to_csv(processed_dir / "gse317012_calibration.csv", index=False)

    # Per-cell export (drop embeddings-heavy columns to keep parquet small)
    obs_cols = [
        c for c in obs.columns if not c.startswith("mk_")
    ]
    obs[obs_cols].to_parquet(
        processed_dir / "gse317012_cells.parquet", index=True
    )

    summary = (
        obs.groupby(["leiden", "cell_type"])
        .size()
        .reset_index(name="n_cells")
        .merge(
            obs.groupby(["leiden", "cell_type"])["viral_response_score"]
            .mean()
            .reset_index(name="mean_viral_response_score"),
            on=["leiden", "cell_type"],
        )
    )
    summary.to_csv(
        processed_dir / "gse317012_cluster_summary.csv", index=False
    )

    model_mapping = build_model_mapping(obs)
    model_mapping.to_csv(
        processed_dir / "gse317012_model_mapping.csv", index=False
    )
    return obs


def build_model_mapping(obs: pd.DataFrame) -> pd.DataFrame:
    """Finding -> ODE state variable / parameter mapping table."""
    epi = obs[obs["infected_call"] != "other"]
    rows = [
        {
            "finding": "Infected-signature epithelial cells show elevated "
            "translation/ribosome module",
            "evidence_column": "mod_translation",
            "model_target": "translation_enhancement (param) / "
            "pathway_translation (state)",
            "measured_shift": float(
                epi.loc[epi["infected_call"] == "infected_signature",
                        "mod_translation"].mean()
                - epi.loc[epi["infected_call"] == "control_epithelial",
                          "mod_translation"].mean()
            ),
            "direction": "up",
        },
        {
            "finding": "Elevated cell-cycle/S-G2M module in infected-signature "
            "epithelial cells",
            "evidence_column": "mod_s_g2m",
            "model_target": "s_phase_bonus, cc_s_phase gate (params) / "
            "CC, DNA (states)",
            "measured_shift": float(
                epi.loc[epi["infected_call"] == "infected_signature",
                        "mod_s_g2m"].mean()
                - epi.loc[epi["infected_call"] == "control_epithelial",
                          "mod_s_g2m"].mean()
            ),
            "direction": "up",
        },
        {
            "finding": "DDR module elevated in infected-signature epithelium",
            "evidence_column": "mod_ddr",
            "model_target": "ddr_enhancement (param) / "
            "pathway_dna_damage_response (state)",
            "measured_shift": float(
                epi.loc[epi["infected_call"] == "infected_signature",
                        "mod_ddr"].mean()
                - epi.loc[epi["infected_call"] == "control_epithelial",
                          "mod_ddr"].mean()
            ),
            "direction": "up",
        },
        {
            "finding": "Interferon/antiviral module elevated",
            "evidence_column": "mod_interferon",
            "model_target": "p_immune_prod, ifn_prod (params) / IFN, AK, "
            "P_immune (states)",
            "measured_shift": float(
                epi.loc[epi["infected_call"] == "infected_signature",
                        "mod_interferon"].mean()
                - epi.loc[epi["infected_call"] == "control_epithelial",
                          "mod_interferon"].mean()
            ),
            "direction": "up",
        },
        {
            "finding": "Antigen presentation module elevated",
            "evidence_column": "mod_antigen_presentation",
            "model_target": "pathway_antigen_presentation (state) / "
            "antigen_presentation_suppression_rate (param)",
            "measured_shift": float(
                epi.loc[epi["infected_call"] == "infected_signature",
                        "mod_antigen_presentation"].mean()
                - epi.loc[epi["infected_call"] == "control_epithelial",
                          "mod_antigen_presentation"].mean()
            ),
            "direction": "up",
        },
        {
            "finding": "Mitochondrial module shift in infected-signature "
            "epithelial cells",
            "evidence_column": "mod_mito",
            "model_target": "mitochondrial_importance (param) / "
            "pathway_mitochondrial_stress (state)",
            "measured_shift": float(
                epi.loc[epi["infected_call"] == "infected_signature",
                        "mod_mito"].mean()
                - epi.loc[epi["infected_call"] == "control_epithelial",
                          "mod_mito"].mean()
            ),
            "direction": "context-dependent",
        },
        {
            "finding": "ECM/wound-healing module shift (epithelial stress)",
            "evidence_column": "mod_ecm_wound",
            "model_target": "pathway_cellular_stress (state)",
            "measured_shift": float(
                epi.loc[epi["infected_call"] == "infected_signature",
                        "mod_ecm_wound"].mean()
                - epi.loc[epi["infected_call"] == "control_epithelial",
                          "mod_ecm_wound"].mean()
            ),
            "direction": "up",
        },
    ]
    return pd.DataFrame(rows)


if __name__ == "__main__":
    raise SystemExit(main())
