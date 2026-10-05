"""Command-line interface for VCM."""

import sys
from pathlib import Path
from typing import Optional

import click
from rich.console import Console
from rich.table import Table

# Windows consoles default to cp1252 and crash on glyphs like "✓"
try:  # pragma: no cover - platform dependent
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except (AttributeError, OSError):
    pass

from vcm.experiments.runner import ExperimentComparator, ExperimentRunner
from vcm.plugins.bacteria.minimal_cell import MinimalCellPlugin
from vcm.plugins.mammalian.cancer_cell import CancerCellPlugin
from vcm.plugins.mammalian.immune_tcell import ImmuneTCellPlugin
from vcm.plugins.virus_host.sars_cov2 import SARSCoV2Plugin
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.plugins.base import PluginRegistry

console = Console()


# Register plugins
def register_default_plugins():
    """Register default plugins."""
    plugins = [
        MinimalCellPlugin(),
        SARSCoV2Plugin(),
        ImmuneTCellPlugin(),
        CancerCellPlugin(),
        BKPolyomavirusPlugin(),
    ]
    for plugin in plugins:
        plugin.register()


@click.group()
@click.version_option(version="0.1.0")
def cli():
    """Virtual Cell Model (VCM) - A modular platform for cellular simulation."""
    register_default_plugins()


@cli.command()
@click.option("--plugin", required=True, help="Plugin ID to use")
@click.option("--config", required=True, type=click.Path(exists=True), help="Path to config file")
@click.option("--output", default=None, help="Output directory (default: from config)")
def run(plugin: str, config: str, output: Optional[str]):
    """Run a simulation experiment."""
    try:
        # Initialize runner
        runner = ExperimentRunner()

        # Load config
        config_path = Path(config)
        experiment_config = runner.load_config(config_path)

        # Override plugin if specified
        if plugin:
            experiment_config.plugin = plugin

        # Override output if specified
        if output:
            experiment_config.output_path = output

        console.print(f"[bold blue]Running experiment:[/bold blue] {experiment_config.experiment_id}")
        console.print(f"Plugin: {experiment_config.plugin}")
        console.print(f"Simulator: {experiment_config.simulator}")
        console.print(f"Simulation length: {experiment_config.simulation_length}")
        console.print(f"Timestep: {experiment_config.timestep}")
        console.print(f"Perturbations: {len(experiment_config.perturbations)}")
        console.print()

        # Run experiment
        with console.status("[bold green]Running simulation..."):
            result = runner.run_experiment(experiment_config)

        console.print("[bold green]✓ Simulation complete![/bold green]")
        console.print(f"Steps simulated: {len(result.steps)}")
        console.print(f"Success: {result.success}")

        # Save results
        output_dir = Path(experiment_config.output_path)
        runner.save_result(result, output_dir)
        console.print(f"[bold green]✓ Results saved to:[/bold green] {output_dir}")

    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.ClickException(str(e))


@cli.command()
@click.option("--plugin", required=True, help="Plugin ID to use")
@click.option("--scenario1", required=True, help="First scenario config file")
@click.option("--scenario2", required=True, help="Second scenario config file")
@click.option("--entity-id", default=None, help="Entity ID to compare trajectories")
@click.option("--entity-type", default="gene", help="Entity type (gene, protein, metabolite, pathway)")
def compare(plugin: str, scenario1: str, scenario2: str, entity_id: Optional[str], entity_type: str):
    """Compare two simulation scenarios."""
    try:
        # Initialize runner and comparator
        runner = ExperimentRunner()
        comparator = ExperimentComparator()

        # Load and run scenario 1
        console.print(f"[bold blue]Running scenario 1:[/bold blue] {scenario1}")
        config1 = runner.load_config(Path(scenario1))
        config1.plugin = plugin
        result1 = runner.run_experiment(config1)
        comparator.add_result("scenario1", result1)
        console.print("[bold green]✓ Scenario 1 complete[/bold green]")

        # Load and run scenario 2
        console.print(f"[bold blue]Running scenario 2:[/bold blue] {scenario2}")
        config2 = runner.load_config(Path(scenario2))
        config2.plugin = plugin
        result2 = runner.run_experiment(config2)
        comparator.add_result("scenario2", result2)
        console.print("[bold green]✓ Scenario 2 complete[/bold green]")
        console.print()

        # Compare final states
        console.print("[bold blue]Comparing final states:[/bold blue]")
        comparison = comparator.compare_states("scenario1", "scenario2")

        # Create table for gene differences
        gene_table = Table(title="Gene Expression Differences")
        gene_table.add_column("Gene ID", style="cyan")
        gene_table.add_column("Difference", style="magenta")

        for gene_id, diff in comparison["gene_differences"].items():
            gene_table.add_row(gene_id, f"{diff:.4f}")

        console.print(gene_table)

        # Create table for protein differences
        protein_table = Table(title="Protein Concentration Differences")
        protein_table.add_column("Protein ID", style="cyan")
        protein_table.add_column("Difference", style="magenta")

        for protein_id, diff in comparison["protein_differences"].items():
            protein_table.add_row(protein_id, f"{diff:.4f}")

        console.print(protein_table)

        # Compare trajectories if entity specified
        if entity_id:
            console.print(f"[bold blue]Comparing trajectories for {entity_type}: {entity_id}[/bold blue]")
            traj_comparison = comparator.compare_trajectories("scenario1", "scenario2", entity_id, entity_type)

            traj_table = Table(title=f"Trajectory Comparison: {entity_id}")
            traj_table.add_column("Timepoint", style="cyan")
            traj_table.add_column("Scenario 1", style="green")
            traj_table.add_column("Scenario 2", style="blue")
            traj_table.add_column("Difference", style="magenta")

            for i, tp in enumerate(traj_comparison["timepoints"]):
                val1 = traj_comparison["trajectory1"][i] if i < len(traj_comparison["trajectory1"]) else "N/A"
                val2 = traj_comparison["trajectory2"][i] if i < len(traj_comparison["trajectory2"]) else "N/A"
                diff = traj_comparison["differences"][i] if i < len(traj_comparison["differences"]) else "N/A"
                traj_table.add_row(f"{tp:.1f}", str(val1), str(val2), str(diff))

            console.print(traj_table)

    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.ClickException(str(e))


@cli.command()
def list_plugins():
    """List all available plugins."""
    register_default_plugins()
    plugins = PluginRegistry.list_plugins()

    console.print("[bold blue]Available Plugins:[/bold blue]")
    table = Table()
    table.add_column("Plugin ID", style="cyan")
    table.add_column("Plugin Name", style="green")
    table.add_column("Description", style="white")

    for plugin_id in plugins:
        plugin = PluginRegistry.get(plugin_id)
        if plugin:
            table.add_row(plugin.plugin_id, plugin.plugin_name, plugin.description)

    console.print(table)


@cli.command()
@click.option("--plugin", required=True, help="Plugin ID to inspect")
def inspect(plugin: str):
    """Inspect a plugin's schema and supported perturbations."""
    register_default_plugins()
    plugin_obj = PluginRegistry.get(plugin)

    if plugin_obj is None:
        console.print(f"[bold red]Plugin '{plugin}' not found[/bold red]")
        raise click.ClickException(f"Plugin '{plugin}' not found")

    console.print(f"[bold blue]Plugin:[/bold blue] {plugin_obj.plugin_name}")
    console.print(f"[bold blue]ID:[/bold blue] {plugin_obj.plugin_id}")
    console.print(f"[bold blue]Description:[/bold blue] {plugin_obj.description}")
    console.print()

    # Get schema
    schema = plugin_obj.get_cell_schema()
    console.print("[bold blue]Cell Schema:[/bold blue]")
    schema_table = Table()
    schema_table.add_column("Property", style="cyan")
    schema_table.add_column("Value", style="green")

    schema_table.add_row("Cell type", schema["cell_type"])
    schema_table.add_row("Number of genes", str(schema["n_genes"]))
    schema_table.add_row("Number of proteins", str(schema["n_proteins"]))
    schema_table.add_row("Number of metabolites", str(schema["n_metabolites"]))
    schema_table.add_row("Number of pathways", str(schema["n_pathways"]))

    console.print(schema_table)
    console.print()

    # Get supported perturbations
    perturbations = plugin_obj.get_supported_perturbations()
    console.print(f"[bold blue]Supported Perturbations ({len(perturbations)}):[/bold blue]")
    pert_table = Table()
    pert_table.add_column("ID", style="cyan")
    pert_table.add_column("Name", style="green")
    pert_table.add_column("Type", style="magenta")
    pert_table.add_column("Target", style="yellow")

    for pert in perturbations:
        pert_table.add_row(pert.id, pert.name, pert.perturbation_type.value, str(pert.target_id))

    console.print(pert_table)


if __name__ == "__main__":
    cli()
