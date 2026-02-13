"""Pipeline runner — orchestrates stage execution in order.

This module provides a reusable :func:`run_pipeline` function that can be
called from the CLI or used programmatically.  It handles:

- Stage ordering and filtering (single stage, from-stage, configured list)
- Input validation before execution
- Skip logic (outputs exist + config hash unchanged)
- State tracking via :class:`~merfish_pipeline.execution.state.PipelineState`
- Metadata writing after each stage
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Any

from merfish_pipeline.config.defaults import VALID_STAGES
from merfish_pipeline.exceptions import StageError
from merfish_pipeline.execution.state import PipelineState
from merfish_pipeline.stages.base import StageResult
from merfish_pipeline.stages.registry import get_stage

logger = logging.getLogger("merfish_pipeline.execution.runner")


def resolve_stages(
    config: Any,
    single_stage: str | None = None,
    from_stage: str | None = None,
) -> list[str]:
    """Determine which stages to run based on CLI options and config.

    Priority:
    1. ``single_stage`` — run exactly one stage.
    2. ``from_stage`` — run from this stage onward within the configured list.
    3. ``config.pipeline.stages`` — the user's configured stage list.

    Raises
    ------
    ValueError
        If no stages can be determined.
    """
    if single_stage:
        return [single_stage]

    if from_stage:
        if from_stage not in VALID_STAGES:
            raise ValueError(f"Unknown stage: {from_stage}")
        idx = VALID_STAGES.index(from_stage)
        configured = config.pipeline.stages or VALID_STAGES
        return [s for s in configured if s in VALID_STAGES and VALID_STAGES.index(s) >= idx]

    stages = config.pipeline.stages
    if not stages:
        raise ValueError(
            "No stages configured.  Use --stage to run a specific stage, "
            "or set pipeline.stages in your experiment config."
        )
    return list(stages)


def run_pipeline(
    config: Any,
    stages: list[str],
    state: PipelineState,
    *,
    force: bool = False,
    dry_run: bool = False,
) -> list[dict[str, Any]]:
    """Execute a list of pipeline stages in order.

    Parameters
    ----------
    config:
        Fully resolved ``PipelineConfig`` instance.
    stages:
        Ordered list of stage names to execute.
    state:
        Pipeline state tracker.
    force:
        If True, re-run stages even when outputs already exist.
    dry_run:
        If True, preview actions without writing data.

    Returns
    -------
    list[dict]
        Per-stage execution report with keys: ``name``, ``status``,
        ``duration_seconds``, ``output_files``, ``error``.

    Raises
    ------
    StageError
        If a stage fails and halts the pipeline.
    """
    reports: list[dict[str, Any]] = []

    for stage_name in stages:
        report: dict[str, Any] = {
            "name": stage_name,
            "status": "pending",
            "duration_seconds": 0.0,
            "output_files": [],
            "error": "",
        }

        # --- Resolve stage class ---
        try:
            stage_cls = get_stage(stage_name)
        except KeyError:
            logger.warning(
                "Stage '%s' is not yet implemented (not registered). Skipping.",
                stage_name,
            )
            report["status"] = "skipped"
            report["error"] = "not registered"
            reports.append(report)
            continue

        stage = stage_cls(config=config, state=state)

        # --- Check per-stage enabled flag ---
        stage_cfg = getattr(config, stage_name, None)
        if stage_cfg is not None and hasattr(stage_cfg, "enabled") and not stage_cfg.enabled:
            logger.info("[%s] Disabled via config (enabled: false), skipping.", stage_name)
            report["status"] = "skipped"
            reports.append(report)
            continue

        # --- Validate inputs ---
        errors = stage.validate_inputs()
        if errors:
            for err in errors:
                logger.error("  [%s] Input error: %s", stage_name, err)
            report["status"] = "failed"
            report["error"] = "; ".join(errors)
            reports.append(report)
            raise StageError(stage_name, f"Input validation failed: {'; '.join(errors)}")

        # --- Check skip ---
        if stage.should_skip(force=force):
            logger.info("[%s] Outputs exist, skipping (use --force to re-run).", stage_name)
            report["status"] = "skipped"
            reports.append(report)
            continue

        # --- Dry run ---
        if dry_run:
            logger.info("[%s] DRY RUN -- would execute stage.", stage_name)
            report["status"] = "dry_run"
            reports.append(report)
            continue

        # --- Execute ---
        logger.info("[%s] Starting...", stage_name)
        state.mark_started(stage_name)
        start_time = datetime.now()

        try:
            result = stage.run(dry_run=False)
        except Exception as exc:
            duration = (datetime.now() - start_time).total_seconds()
            state.mark_failed(stage_name, str(exc))
            logger.error("[%s] FAILED after %.1fs: %s", stage_name, duration, exc)
            report["status"] = "failed"
            report["error"] = str(exc)
            report["duration_seconds"] = round(duration, 3)
            reports.append(report)
            raise StageError(stage_name, str(exc), original=exc) from exc

        duration = (datetime.now() - start_time).total_seconds()

        if result.status == "failed":
            state.mark_failed(stage_name, result.error)
            logger.error("[%s] FAILED: %s", stage_name, result.error)
            report["status"] = "failed"
            report["error"] = result.error
            report["duration_seconds"] = round(duration, 3)
            reports.append(report)
            raise StageError(stage_name, result.error)

        state.mark_completed(stage_name, result)
        logger.info(
            "[%s] Completed in %.1fs. Outputs: %d files.",
            stage_name,
            duration,
            len(result.output_files),
        )

        report["status"] = result.status
        report["output_files"] = result.output_files
        report["duration_seconds"] = round(duration, 3)
        reports.append(report)

    return reports
