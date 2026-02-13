"""``merfish-pipe status`` — show pipeline execution status."""

from __future__ import annotations

from pathlib import Path

import click

from merfish_pipeline.execution.state import PipelineState


@click.command()
@click.argument("experiment", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--profile",
    default="local",
    show_default=True,
    help="Execution profile (to resolve output_dir).",
)
def status(experiment: Path, profile: str) -> None:
    """Show stage completion status for an experiment.

    EXPERIMENT is the path to the experiment YAML config file.
    Reads the pipeline state file from the configured output directory.
    """
    from merfish_pipeline.config.loader import load_pipeline_config

    try:
        cfg = load_pipeline_config(experiment, profile=profile)
    except Exception as exc:
        click.echo(f"Config error: {exc}", err=True)
        raise SystemExit(1)

    state = PipelineState.load_or_create(cfg.paths.output_dir)
    summary = state.summary()

    click.echo(f"Experiment: {cfg.experiment.name}")
    click.echo(f"Output dir: {cfg.paths.output_dir}")
    click.echo()

    stages_status = summary.get("stages", {})
    if not stages_status:
        click.echo("No stages have been run yet.")
        return

    # Column widths
    max_name = max(len(n) for n in stages_status)
    click.echo(f"{'Stage':<{max_name + 2}} Status")
    click.echo("-" * (max_name + 14))

    status_symbols = {
        "completed": "done",
        "running": "running",
        "failed": "FAILED",
    }

    for name, st in stages_status.items():
        display = status_symbols.get(st, st)
        click.echo(f"  {name:<{max_name}}  {display}")

    click.echo()
    counts = summary.get("counts", {})
    parts = [f"{v} {k}" for k, v in counts.items()]
    click.echo(f"Summary: {', '.join(parts)}")
