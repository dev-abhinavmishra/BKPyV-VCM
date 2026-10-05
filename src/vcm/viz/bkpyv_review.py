"""Publication-style BKPyV output bundle.

The bundle is deliberately built from the ODE result rather than hand-drawn
clinical curves. It produces figures, a long-form trajectory table, scenario
summaries, and a machine-readable review manifest so another researcher can
inspect exactly what was plotted and what remains hypothetical.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vcm.clinical.viral_load_mapper import ViralLoadMapper
from vcm.core.models import Perturbation, PerturbationType
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator

PALETTE = {
    "archetype": "#2B6CB0",
    "rearranged": "#C05621",
    "viral": "#7B2CBF",
    "host": "#2F855A",
    "immune": "#4A5568",
    "reference": "#718096",
    "screening": "#D69E2E",
    "treatment": "#C53030",
    "paper": "#1A202C",
}

# Chart contract: the bundle contains a temporal trajectory (trend), two
# small-multiple comparisons, and a parameter-by-value sensitivity matrix.
# All figures are generated from the same ODE result tables and carry a
# neutral title plus a caveat-bearing subtitle.
REVIEW_SIMULATION_CONFIG = {
    "ode_solver": "LSODA",
    "rtol": 1e-6,
    "atol": 1e-8,
    "max_step": 1.0,
    "timestep": 1.0,
}


def _scenario_result(nccr_variant: str, drug: str = "none", days: float = 60.0, timestep: float = 1.0):
    plugin = BKPolyomavirusPlugin()
    initial_state = plugin.create_initial_state({"nccr_variant": nccr_variant})
    perturbations = [
        Perturbation(
            id="bkpyv_infection",
            name="BKPyV infection",
            perturbation_type=PerturbationType.VIRAL_INFECTION,
            magnitude=1.0,
            timing=0.0,
        )
    ]
    if drug == "tacrolimus":
        perturbations.append(Perturbation(id="tacrolimus", name="Tacrolimus", perturbation_type=PerturbationType.DRUG_TREATMENT, target_id="FKBP1A", magnitude=1.0, timing=0.0))
    elif drug == "sirolimus":
        perturbations.append(Perturbation(id="sirolimus", name="Sirolimus", perturbation_type=PerturbationType.DRUG_TREATMENT, target_id="MTOR", magnitude=1.0, timing=0.0))
    simulator_config = dict(REVIEW_SIMULATION_CONFIG)
    simulator_config["max_step"] = timestep
    simulator_config["timestep"] = timestep
    simulator = BKPyVODESimulator(simulator_config)
    return simulator.simulate(initial_state, perturbations=perturbations, n_steps=int(days / timestep), timestep=timestep)


def _trajectory_frame(result, scenario: str) -> pd.DataFrame:
    """Convert simulator steps to a stable, reviewable table schema."""
    columns = [
        "time_days", "scenario", "nccr_variant", "virtual_viral_load",
        "plasma_copies_per_ml", "large_T_antigen", "host_DNA_synthesis",
        "cell_cycle_position", "intracellular_replication_flux",
        "viral_production_rate", "immune_control_index",
        "vp1_capsid_expression", "infection_status",
    ]
    rows = []
    mapper = ViralLoadMapper()
    for step in result.steps:
        state = step.cell_state
        meta = getattr(state, "metadata", {}) or {}
        pathway = meta.get("pathway_activities", {}) or {}
        genes = getattr(state, "genes", {}) or {}
        lt_gene = genes.get("viral_LT")
        vp1_gene = genes.get("viral_VP1")
        rows.append({
            "time_days": float(step.timestamp),
            "scenario": scenario,
            "nccr_variant": meta.get("nccr_variant", "archetype"),
            "virtual_viral_load": float(meta.get("viral_load", 0.0)),
            "plasma_copies_per_ml": float(mapper.normalized_to_copies(meta.get("viral_load", 0.0))),
            # ``t_antigen_level`` is a categorical label in the simulator
            # metadata; the gene expression value is the numeric trace.
            "large_T_antigen": float(getattr(lt_gene, "expression_level", 0.0)),
            "host_DNA_synthesis": float(pathway.get("dna_replication", 0.0)),
            "cell_cycle_position": float(pathway.get("cell_cycle", 0.0)),
            "intracellular_replication_flux": float(meta.get("intracellular_replication_flux", 0.0)),
            "viral_production_rate": float(meta.get("viral_production_rate", 0.0)),
            "immune_control_index": float(meta.get("immune_control_index", 0.0)),
            "vp1_capsid_expression": float(getattr(vp1_gene, "expression_level", 0.0)),
            "infection_status": meta.get("infection_status", "uninfected"),
        })
    return pd.DataFrame(rows, columns=columns)


def _style(ax, title: str, subtitle: str, ylabel: str) -> None:
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", color=PALETTE["paper"], pad=18)
    ax.text(0, 1.01, subtitle, transform=ax.transAxes, fontsize=8.5, color=PALETTE["reference"], va="bottom")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#E2E8F0", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#CBD5E0")


def plot_mechanistic_trajectory(df: pd.DataFrame, output_path: Path) -> None:
    """Plot one scenario with biology and clinical bridge separated."""
    if df.empty:
        raise ValueError("Cannot plot a mechanistic trajectory from an empty result")
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True, constrained_layout=True)
    x = df["time_days"]
    axes[0].plot(x, df["plasma_copies_per_ml"].clip(lower=1), color=PALETTE["viral"], lw=2.2)
    axes[0].axhline(1_000, color=PALETTE["screening"], ls="--", lw=1.2, label="Screening reference")
    axes[0].axhline(10_000, color=PALETTE["treatment"], ls="--", lw=1.2, label="Treatment reference")
    axes[0].set_yscale("log")
    _style(axes[0], "Clinical bridge: simulated plasma BKPyV signal", "Mapper output shown for orientation; not a patient-calibrated prediction", "copies/mL")
    axes[0].legend(frameon=False, ncol=2, loc="upper left")

    axes[1].plot(x, df["large_T_antigen"], color=PALETTE["viral"], lw=2, label="Large T antigen")
    axes[1].plot(x, df["host_DNA_synthesis"], color=PALETTE["host"], lw=2, label="Host DNA synthesis")
    axes[1].plot(x, df["cell_cycle_position"], color=PALETTE["immune"], lw=1.8, ls="--", label="Cell-cycle position")
    _style(axes[1], "Mechanistic gate: host state precedes LT accumulation", "State variables are normalized/model units; inspect the gate before interpreting output", "model units")
    axes[1].legend(frameon=False, ncol=3, loc="upper left")

    axes[2].plot(x, df["intracellular_replication_flux"], color=PALETTE["viral"], lw=2, label="Genome-copying flux")
    axes[2].plot(x, df["viral_production_rate"], color=PALETTE["rearranged"], lw=2, label="Production rate")
    axes[2].plot(x, df["immune_control_index"], color=PALETTE["immune"], lw=1.8, ls=":", label="Immune-control index")
    _style(axes[2], "Mechanism separation: replication, production, and immune control", "Tacrolimus should move immune control first; it is not a direct genome-copying multiplier", "model units")
    axes[2].set_xlabel("Time (days)")
    axes[2].legend(frameon=False, ncol=3, loc="upper left")
    fig.suptitle("BKPyV virtual-cell trajectory review", x=0.08, ha="left", fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_nccr_comparison(frames: Mapping[str, pd.DataFrame], output_path: Path) -> None:
    """Compare archetype/rearranged NCCR as a falsifiable hypothesis."""
    if not frames or any(frame.empty for frame in frames.values()):
        raise ValueError("NCCR comparison requires non-empty trajectories")
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    metrics = [
        ("large_T_antigen", "Large T antigen", "model units"),
        ("vp1_capsid_expression", "VP1 capsid expression", "model units"),
        ("intracellular_replication_flux", "Genome-copying flux", "model units"),
        ("plasma_copies_per_ml", "Clinical bridge", "copies/mL"),
    ]
    for ax, (column, title, ylabel) in zip(axes.flat, metrics):
        for variant, frame in frames.items():
            y = frame[column].clip(lower=1e-6) if column == "plasma_copies_per_ml" else frame[column]
            ax.plot(frame["time_days"], y, lw=2, color=PALETTE[variant], label=variant.capitalize())
        if column == "plasma_copies_per_ml":
            ax.set_yscale("log")
        _style(ax, title, "Same host/drug context; NCCR is the varied factor", ylabel)
        ax.set_xlabel("Time (days)")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("NCCR scenario comparison", x=0.08, ha="left", fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_drug_mechanisms(frames: Mapping[str, pd.DataFrame], output_path: Path) -> None:
    """Show drug effects by mechanism, not as a single replication multiplier."""
    if not frames or any(frame.empty for frame in frames.values()):
        raise ValueError("Drug comparison requires non-empty trajectories")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), sharex=True, constrained_layout=True)
    metrics = [
        ("immune_control_index", "Immune-control index", "higher = more control"),
        ("intracellular_replication_flux", "Genome-copying flux", "intracellular model units"),
        ("viral_production_rate", "System-level production rate", "model units"),
    ]
    colors = {"none": PALETTE["reference"], "tacrolimus": PALETTE["rearranged"], "sirolimus": PALETTE["host"]}
    for ax, (column, title, ylabel) in zip(axes, metrics):
        for drug, frame in frames.items():
            ax.plot(frame["time_days"], frame[column], lw=2, color=colors[drug], label=drug.capitalize())
        _style(ax, title, "Archetype NCCR; drug effect is pathway-specific", ylabel)
        ax.set_xlabel("Time (days)")
    axes[0].legend(frameon=False)
    fig.suptitle("Drug mechanism comparison", x=0.08, ha="left", fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _sensitivity_table(nccr_variant: str, output_dir: Path, days: float) -> pd.DataFrame:
    """One-at-a-time sensitivity on interpretable review outcomes."""
    parameters = {
        "t_antigen_replication_threshold": [0.3, 0.4, 0.5, 0.6, 0.7],
        "dna_replication_coupling": [0.5, 0.65, 0.8, 0.9, 1.0],
        "nccr_early_expression_multiplier": [0.6, 0.8, 1.0, 1.4, 2.0],
        "nccr_capsid_expression_multiplier": [0.3, 0.5, 0.75, 1.0, 1.2],
    }
    rows = []
    plugin = BKPolyomavirusPlugin()
    for parameter, values in parameters.items():
        for value in values:
            state = plugin.create_initial_state({"nccr_variant": nccr_variant})
            sim = BKPyVODESimulator({"ode_solver": "LSODA", parameter: value})
            infection = Perturbation(id="infection", name="BKPyV infection", perturbation_type=PerturbationType.VIRAL_INFECTION, magnitude=1.0, timing=0.0)
            result = sim.simulate(state, perturbations=[infection], n_steps=int(days), timestep=1.0)
            frame = _trajectory_frame(result, parameter)
            positive = frame.loc[frame["large_T_antigen"] >= 0.5, "time_days"]
            rows.append({
                "parameter": parameter,
                "value": value,
                "peak_production_rate": frame["viral_production_rate"].max(),
                "peak_replication_flux": frame["intracellular_replication_flux"].max(),
                "time_to_LT_threshold_days": float(positive.iloc[0]) if len(positive) else np.nan,
            })
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "sensitivity_results.csv", index=False)
    return table


def plot_sensitivity_heatmap(table: pd.DataFrame, output_path: Path) -> None:
    if table.empty:
        raise ValueError("Sensitivity heatmap requires at least one result row")
    pivot = table.pivot(index="parameter", columns="value", values="peak_production_rate")
    fig, ax = plt.subplots(figsize=(11, 4.8), constrained_layout=True)
    image = ax.imshow(pivot.values, aspect="auto", cmap="Purples")
    ax.set_xticks(range(len(pivot.columns)), [f"{v:g}" for v in pivot.columns])
    ax.set_yticks(range(len(pivot.index)), [p.replace("_", " ").title() for p in pivot.index])
    ax.set_xlabel("Parameter value")
    ax.set_title("Parameter sensitivity: peak viral production rate", loc="left", fontsize=13, fontweight="bold", pad=20)
    ax.text(0, 1.01, "One-at-a-time sweep; darker cells indicate larger model output, not greater evidence", transform=ax.transAxes, fontsize=8.5, color=PALETTE["reference"], va="bottom")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.iloc[i, j]
            ax.text(j, i, f"{value:.2g}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, label="Peak production rate (model units)")
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_review_bundle(output_dir: str | Path = "outputs/bkpyv_review", days: float = 60.0) -> Dict[str, str]:
    """Generate all reviewer-facing outputs and return their paths."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    archetype = _trajectory_frame(_scenario_result("archetype", days=days), "archetype")
    rearranged = _trajectory_frame(_scenario_result("rearranged", days=days), "rearranged")
    combined = pd.concat([archetype, rearranged], ignore_index=True)
    combined.to_csv(output_path / "trajectory_data.csv", index=False)
    plot_mechanistic_trajectory(archetype, output_path / "mechanistic_trajectory.png")
    plot_nccr_comparison({"archetype": archetype, "rearranged": rearranged}, output_path / "nccr_comparison.png")
    drug_frames = {
        drug: _trajectory_frame(_scenario_result("archetype", drug=drug, days=days), drug)
        for drug in ("none", "tacrolimus", "sirolimus")
    }
    pd.concat(drug_frames.values(), ignore_index=True).to_csv(output_path / "drug_trajectories.csv", index=False)
    plot_drug_mechanisms(drug_frames, output_path / "drug_mechanisms.png")
    sensitivity = _sensitivity_table("archetype", output_path, days=min(days, 30.0))
    plot_sensitivity_heatmap(sensitivity, output_path / "sensitivity_heatmap.png")
    summary = combined.groupby("scenario").agg(
        peak_plasma_copies_per_ml=("plasma_copies_per_ml", "max"),
        peak_replication_flux=("intracellular_replication_flux", "max"),
        peak_production_rate=("viral_production_rate", "max"),
        peak_LT=("large_T_antigen", "max"),
        peak_VP1=("vp1_capsid_expression", "max"),
    ).reset_index()
    summary.to_csv(output_path / "scenario_summary.csv", index=False)
    manifest = {
        "title": "BKPyV virtual-cell review bundle",
        "model_version": "S-phase-gated ODE with NCCR scenario axis",
        "simulation_days": days,
        "days": days,
        "timestep": REVIEW_SIMULATION_CONFIG["timestep"],
        "ode_solver": REVIEW_SIMULATION_CONFIG["ode_solver"],
        "rtol": REVIEW_SIMULATION_CONFIG["rtol"],
        "atol": REVIEW_SIMULATION_CONFIG["atol"],
        "max_step": REVIEW_SIMULATION_CONFIG["max_step"],
        "nccr_variant": ["archetype", "rearranged"],
        "scenarios": ["archetype", "rearranged"],
        "figures": ["mechanistic_trajectory.png", "nccr_comparison.png", "drug_mechanisms.png", "sensitivity_heatmap.png"],
        "tables": ["trajectory_data.csv", "drug_trajectories.csv", "scenario_summary.csv", "sensitivity_results.csv"],
        "interpretation": "Hypothesis-generating; clinical bridge is not patient calibration.",
        "interpretation_caveats": [
            "Plasma copies/mL are a population mapper from normalized model output, not patient calibration.",
            "NCCR archetype/rearranged outputs are falsifiable scenario comparisons, not prevalence estimates.",
            "Drug curves separate immune-control and S-phase mechanisms; they are not treatment recommendations.",
            "Sensitivity results are one-at-a-time parameter sweeps and do not represent posterior uncertainty.",
        ],
    }
    (output_path / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return {key: str(output_path / filename) for key, filename in {"trajectory": "mechanistic_trajectory.png", "nccr": "nccr_comparison.png", "drugs": "drug_mechanisms.png", "sensitivity": "sensitivity_heatmap.png", "trajectory_data": "trajectory_data.csv", "summary": "scenario_summary.csv", "manifest": "manifest.json"}.items()}
