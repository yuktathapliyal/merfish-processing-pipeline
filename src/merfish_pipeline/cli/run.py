"""``merfish-pipe run`` -- execute pipeline stages."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from merfish_pipeline.config.defaults import VALID_STAGES


@click.command()
@click.argument("experiment", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--profile",
    default="local",
    show_default=True,
    help="Execution profile name (e.g. 'local', 'slurm').",
)
@click.option(
    "--stage",
    "single_stage",
    default=None,
    type=click.Choice(VALID_STAGES, case_sensitive=True),
    help="Run a single stage instead of the configured stage list.",
)
@click.option(
    "--from-stage",
    default=None,
    type=click.Choice(VALID_STAGES, case_sensitive=True),
    help="Resume from this stage onward (skips earlier stages).",
)
@click.option("--dry-run", is_flag=True, help="Preview what would be done without executing.")
@click.option("--force", is_flag=True, help="Re-run stages even if outputs already exist.")
@click.option("--workers", type=int, default=None, help="Override max_workers.")
@click.option(
    "--slurm", "use_slurm", is_flag=True, help="Override execution mode to SLURM."
)
@click.option("-v", "--verbose", is_flag=True, help="Increase log verbosity to DEBUG.")
@click.option(
    "--slurm-worker",
    is_flag=True,
    hidden=True,
    help="Internal flag set by generated SLURM scripts to execute stages directly.",
)
def run(
    experiment: Path,
    profile: str,
    single_stage: str | None,
    from_stage: str | None,
    dry_run: bool,
    force: bool,
    workers: int | None,
    use_slurm: bool,
    verbose: bool,
    slurm_worker: bool,
) -> None:
    """Run pipeline stages for an experiment.

    EXPERIMENT is the path to the experiment YAML config file.
    """
    from merfish_pipeline.cli.utils import StatusBanner

    # --- Build CLI overrides ---
    overrides: dict = {}
    if workers is not None:
        overrides.setdefault("execution", {})["max_workers"] = workers
    if use_slurm:
        overrides.setdefault("execution", {})["mode"] = "slurm"
        if profile == "local":
            profile = "slurm"
    if dry_run:
        overrides.setdefault("pipeline", {})["dry_run"] = True
    if force:
        overrides.setdefault("pipeline", {})["force"] = True

    # --- Load config and heavy dependencies (with spinner for feedback) ---
    with StatusBanner("Loading pipeline"):
        from merfish_pipeline.config.loader import load_pipeline_config
        from merfish_pipeline.exceptions import StageError
        from merfish_pipeline.execution.runner import resolve_stages, run_pipeline
        from merfish_pipeline.execution.slurm import generate_slurm_script
        from merfish_pipeline.execution.state import PipelineState
        from merfish_pipeline.logging_config import setup_logging

        # Ensure all stages are registered.
        import merfish_pipeline.stages  # noqa: F401

        try:
            config = load_pipeline_config(
                experiment, profile=profile, overrides=overrides or None
            )
        except Exception as exc:
            raise click.ClickException(f"Config error: {exc}") from exc

    # --- Setup logging ---
    log_level = "DEBUG" if verbose else "INFO"
    logger = setup_logging(
        log_dir=config.logs_dir,
        experiment_name=config.experiment.name,
        level=log_level,
    )

    logger.info(
        "merfish-pipe run: experiment=%s microscope=%s profile=%s",
        config.experiment.name,
        config.experiment.microscope,
        profile,
    )

    # --- Determine which stages to run ---
    try:
        stages_to_run = resolve_stages(config, single_stage=single_stage, from_stage=from_stage)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    logger.info("Stages to run: %s", ", ".join(stages_to_run))

    # --- SLURM mode: generate submission script and exit ---
    # When --slurm-worker is set (by generated SLURM jobs), skip script
    # generation and fall through to local execution of the single stage.
    if config.execution.mode == "slurm" and not slurm_worker:
        script_path = generate_slurm_script(
            config=config,
            stages=stages_to_run,
            experiment_yaml=experiment,
        )
        click.echo(f"SLURM submission script generated: {script_path}")
        click.echo(f"Submit with: bash {script_path}")
        return

    # --- Local mode: execute stages ---
    state = PipelineState.load_or_create(config.paths.output_dir)
    state.experiment = config.experiment.name

    try:
        reports = run_pipeline(
            config=config,
            stages=stages_to_run,
            state=state,
            force=config.pipeline.force,
            dry_run=config.pipeline.dry_run,
        )
    except StageError as exc:
        logger.error("Pipeline stopped: %s", exc)
        sys.exit(1)

    # --- Summary ---
    completed = sum(1 for r in reports if r["status"] == "completed")
    skipped = sum(1 for r in reports if r["status"] == "skipped")
    failed = sum(1 for r in reports if r["status"] == "failed")

    logger.info(
        "Pipeline finished: %d completed, %d skipped, %d failed.",
        completed,
        skipped,
        failed,
    )

    if failed:
        sys.exit(1)
