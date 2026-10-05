#!/usr/bin/env python3
"""Compare the BKPyV ODE model against longitudinal patient viral-load series.

WHY THIS EXISTS
---------------
A collaborating group has offered longitudinal first-year post-transplant
BKPyV series (urine AND plasma, from their measured patient data) to compare
against model scenarios. No patient data ships with this repository: the script
runs only when a real CSV is supplied and simply declines to guess. This keeps
any future claim of "agreement with clinical data" honest: it has to come from
an actual file.

EXPECTED INPUT FORMAT (CSV), one row per sample
-----------------------------------------------
    patient_id          int or str (any cohort label)
    day_post_transplant float       days since transplant
    plasma_copies_ml    float       plasma BKPyV DNAemia (0 if none measured)
    urine_copies_ml     float       optional; urine BKPyV DNA (may be blank)
    drug_group          str         optional; 'tacrolimus' | 'sirolimus' | other

Example template: run with `--make-template` to write
``data/research/patient_series_template.csv``.

WHAT IT REPORTS
---------------
Per patient + per group: first day above 1,000 and 10,000 cp/mL, peak,
day of peak, post-peak log10 decay/week. Model scenarios (infection,
tacrolimus, sirolimus) are run through the ODE engine and scored with the same
metrics. The output juxtaposes model vs data side by side; it does NOT compute
a fitted "match" and should not be summarized that way.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcm.clinical.viral_load_mapper import ViralLoadMapper  # noqa: E402

TEMPLATE_HEADER = (
    "patient_id,day_post_transplant,plasma_copies_ml,urine_copies_ml,drug_group\n"
)


def make_template(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(TEMPLATE_HEADER)
    print(f"Wrote empty template: {path}")


def series_metrics(days: np.ndarray, values: np.ndarray) -> dict:
    """Mechanical metrics for one load series (days, copies/mL)."""
    order = np.argsort(days)
    days, values = days[order], np.maximum(values[order], 0.0)
    above_1k = days[values >= 1_000]
    above_10k = days[values >= 10_000]
    peak_i = int(np.argmax(values))
    peak = float(values[peak_i])
    peak_day = float(days[peak_i])
    # Post-peak log decay per week (only where > 0)
    tail_d = days[peak_i:]
    tail_v = values[peak_i:]
    positive = tail_v > 0
    decay = None
    if positive.sum() >= 2 and tail_v[positive][0] > tail_v[positive][-1]:
        dto = tail_d[positive]
        lv = np.log10(tail_v[positive])
        sinks = np.polyfit(dto, lv, 1)
        decay = float(-sinks[0] * 7.0)
    return {
        "peak_copies_ml": peak,
        "day_of_peak": peak_day,
        "day_first_ge_1000": (float(above_1k[0]) if len(above_1k) else None),
        "day_first_ge_10000": (float(above_10k[0]) if len(above_10k) else None),
        "post_peak_log10_decay_per_week": decay,
        "n_samples": int(len(days)),
    }


def model_reference() -> dict:
    """Same metrics for the model's canonical situations (52-week horizon)."""
    mapper = ViralLoadMapper()
    out = {}
    for scenario in ("infection", "tacrolimus", "sirolimus"):
        df = mapper.simulate_clinical_trajectory(scenario, weeks=52)
        m = series_metrics(
            (df["week"] * 7.0).to_numpy(), df["copies_per_ml"].to_numpy()
        )
        m["bridge"] = "assumption-labelled anchors (not patient-calibrated)"
        out[scenario] = m
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--input",
        type=Path,
        default=None,
        help="CSV of longitudinal patient series (see module docstring)",
    )
    ap.add_argument(
        "--make-template",
        type=Path,
        nargs="?",
        const=Path("data/research/patient_series_template.csv"),
        help="Write an empty input CSV template and exit",
    )
    ap.add_argument(
        "--output-dir", type=Path, default=Path("outputs/patient_compare")
    )
    args = ap.parse_args()

    if args.make_template is not None:
        make_template(args.make_template)
        return 0

    args.output_dir.mkdir(parents=True, exist_ok=True)

    model = model_reference()
    (args.output_dir / "model_reference_metrics.json").write_text(
        json.dumps(model, indent=2)
    )
    print("Model reference metrics written (infection/tacrolimus/sirolimus).")

    if args.input is None or not args.input.exists():
        print(
            "No patient CSV supplied (use --input). Model reference only; "
            "patient comparison requires real data. Exiting gracefully."
        )
        return 0

    df = pd.read_csv(args.input)
    required = {"patient_id", "day_post_transplant", "plasma_copies_ml"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing columns: {sorted(missing)}")

    rows = []
    for pid, grp in df.groupby("patient_id"):
        m = series_metrics(
            grp["day_post_transplant"].to_numpy(dtype=float),
            grp["plasma_copies_ml"].to_numpy(dtype=float),
        )
        m["patient_id"] = str(pid)
        m["drug_group"] = str(grp["drug_group"].iloc[0]) if "drug_group" in grp else "unknown"
        rows.append(m)

    patient_df = pd.DataFrame(rows)
    patient_df.to_csv(args.output_dir / "patient_metrics.csv", index=False)

    # Side-by-side report blocks keyed by drug group when available
    report = {"model_reference": model, "patient_metrics": rows}
    (args.output_dir / "comparison_report.json").write_text(json.dumps(report, indent=2))
    print(f"Patient metrics: {args.output_dir / 'patient_metrics.csv'}")
    print(f"Comparison report: {args.output_dir / 'comparison_report.json'}")
    print(
        "Reminder: this report juxtaposes model outputs with patient series; "
        "it is a first-pass comparison, NOT model validation."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
