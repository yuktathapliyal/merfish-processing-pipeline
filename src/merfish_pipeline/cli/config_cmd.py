"""``merfish-pipe config`` — configuration utilities."""

from __future__ import annotations

import json
from pathlib import Path

import click
import yaml

from merfish_pipeline.config.defaults import VALID_MICROSCOPES


@click.group()
def config() -> None:
    """Configuration helpers: init, validate, show."""


@config.command()
@click.option(
    "--microscope",
    required=True,
    type=click.Choice(VALID_MICROSCOPES, case_sensitive=False),
    help="Microscope type to generate the template for.",
)
@click.option(
    "--output",
    "-o",
    default=None,
    type=click.Path(path_type=Path),
    help="Output file path. Defaults to stdout.",
)
@click.option(
    "--detailed",
    is_flag=True,
    default=False,
    help="Generate a richly-commented template with full documentation for every field.",
)
def init(microscope: str, output: Path | None, detailed: bool) -> None:
    """Generate an experiment config template.

    Produces a YAML template pre-filled with the selected microscope's defaults.
    Edit the paths and stage list to match your experiment.

    Use --detailed for a version with extensive comments explaining every field
    and stage. Without --detailed, you get a minimal template that's quick to
    fill in.
    """
    if detailed:
        content = _load_detailed_template(microscope.lower())
    else:
        template = _build_template(microscope.lower())

        yaml_text = yaml.dump(template, default_flow_style=False, sort_keys=False)

        # Add a header comment
        header = (
            f"# merFISH Pipeline — Experiment Configuration Template\n"
            f"# Microscope: {microscope.upper()}\n"
            f"#\n"
            f"# Edit the paths and settings below for your experiment.\n"
            f"# Run: merfish-pipe config validate <this_file>\n"
            f"# Tip: use --detailed for a fully-commented version.\n"
            f"#\n\n"
        )

        content = header + yaml_text

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(content, encoding="utf-8")
        click.echo(f"Template written to {output}")
    else:
        click.echo(content)


@config.command()
@click.argument("experiment", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--profile",
    default="local",
    show_default=True,
    help="Execution profile to merge.",
)
def validate(experiment: Path, profile: str) -> None:
    """Validate and merge experiment config without running anything.

    Reports errors in the config and exits with code 1 on failure.
    """
    from merfish_pipeline.config.loader import load_pipeline_config

    try:
        cfg = load_pipeline_config(experiment, profile=profile)
    except Exception as exc:
        click.echo(f"INVALID: {exc}", err=True)
        raise SystemExit(1)

    click.echo(f"VALID: experiment={cfg.experiment.name} microscope={cfg.experiment.microscope}")
    click.echo(f"  output_dir: {cfg.paths.output_dir}")
    click.echo(f"  raw_data_dir: {cfg.paths.raw_data_dir}")
    click.echo(f"  execution: {cfg.execution.mode} (workers={cfg.execution.max_workers})")
    if cfg.pipeline.stages:
        click.echo(f"  stages: {', '.join(cfg.pipeline.stages)}")
    else:
        click.echo("  stages: (none configured)")


@config.command()
@click.argument("experiment", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--profile",
    default="local",
    show_default=True,
    help="Execution profile to merge.",
)
@click.option("--json", "as_json", is_flag=True, help="Output as JSON instead of YAML.")
def show(experiment: Path, profile: str, as_json: bool) -> None:
    """Show the fully resolved (merged) config.

    Loads all three layers, merges them, and prints the result.
    """
    from merfish_pipeline.config.loader import load_pipeline_config

    try:
        cfg = load_pipeline_config(experiment, profile=profile)
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        raise SystemExit(1)

    data = cfg.model_dump(mode="json")

    if as_json:
        click.echo(json.dumps(data, indent=2, default=str))
    else:
        click.echo(yaml.dump(data, default_flow_style=False, sort_keys=False))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_detailed_template(microscope: str) -> str:
    """Load the richly-commented XP_template.yaml and patch the microscope."""
    # Locate configs/experiments/ relative to the package root
    pkg_root = Path(__file__).resolve().parents[2]  # src/merfish_pipeline/cli -> src
    template_path = pkg_root.parent / "configs" / "experiments" / "XP_template.yaml"

    if not template_path.exists():
        raise click.ClickException(
            f"Detailed template not found at {template_path}. "
            f"Use config init without --detailed, or ensure the configs/ "
            f"directory is present in the repository."
        )

    text = template_path.read_text(encoding="utf-8")
    # Replace the microscope placeholder in the template
    text = text.replace('microscope: "oni"', f'microscope: "{microscope}"')
    text = text.replace("microscope: oni", f'microscope: "{microscope}"')
    return text


def _build_template(microscope: str) -> dict:
    """Build a template dict for the given microscope."""
    return {
        "experiment": {
            "name": "EXPERIMENT_NAME",
            "microscope": microscope,
        },
        "paths": {
            "raw_data_dir": "/path/to/raw/data",
            "output_dir": "/path/to/output",
        },
        "raw_data": {
            "bead_channel_folder": (
                "488nm, Raw" if microscope != "andor" else "488nm"
            ),
            "data_org_template": None,
            "stage_file": None,
        },
        "merlin": {
            "codebook_template": "/path/to/codebook.csv",
            "analysis_template": None,
            "cores": 100,
        },
        "focus_qc": {
            "sigma": 1.0,
            "ksize": 3,
        },
        "stitch": {
            "group_by": "ir",
        },
        "inspect_positions": {
            "log_file": None,
        },
        "reregistration": {
            "enabled": False,
        },
        "filter_barcodes": {
            "enabled": False,
            "mode": "any",
        },
        "correlation": {
            "enabled": False,
            "bulk_file": None,
        },
        "optimize_correlation": {
            "enabled": False,
            "correlation_threshold": 0.45,
            "n_attempts": 5,
        },
        "segmentation": {
            "enabled": False,
            "aligned_images_dir": None,
            "mode": "3d",
            "nuclei_bit": 0,
            "total_bits": 16,
            "exclude_bits": [],
            "model_type": "cpsam",
            "diameter": None,
            "median_kernel": 3,
            "batch_size": 8,
            "stitch_threshold": 0.5,
            "flow_threshold": 0.4,
            "cellprob_threshold": 0.0,
            "reference_z_slice": None,
            "z_indexing": 1,
        },
        "cell_assignment": {
            "enabled": False,
            "masks_dir": None,
            "crop_margin": 0,
        },
        "barcode_qc": {
            "enabled": False,
            "top_n_genes": 20,
            "spatial_plots_enabled": True,
            "spatial_plots_columns": 3,
        },
        "anndata_export": {
            "enabled": False,
            "min_barcodes_per_cell": 0,
            "exclude_blanks": True,
        },
        "spatial_visualization": {
            "enabled": False,
            "marker_size": 2,
            "max_points": None,
        },
        "pipeline": {
            "stages": ["index", "stitch", "focus_qc", "inspect_positions"],
        },
    }
