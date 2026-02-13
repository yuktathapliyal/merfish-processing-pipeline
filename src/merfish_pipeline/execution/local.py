"""Local execution backend.

Provides a thin wrapper around :func:`~merfish_pipeline.execution.runner.run_pipeline`
for local (non-SLURM) execution.  In local mode stages run sequentially
in the current process; parallelism happens *within* each stage via
thread/process pools controlled by ``config.execution.max_workers``.
"""

from __future__ import annotations

import logging
from typing import Any

from merfish_pipeline.execution.runner import run_pipeline
from merfish_pipeline.execution.state import PipelineState

logger = logging.getLogger("merfish_pipeline.execution.local")


def run_local(
    config: Any,
    stages: list[str],
    state: PipelineState,
    *,
    force: bool = False,
    dry_run: bool = False,
) -> list[dict[str, Any]]:
    """Execute pipeline stages locally (sequential, in-process).

    This is the default execution mode.  Each stage runs sequentially
    in the current Python process.  Intra-stage parallelism (e.g.
    parallel TIFF I/O) uses the worker count from
    ``config.execution.max_workers``.

    Parameters and return value are identical to
    :func:`~merfish_pipeline.execution.runner.run_pipeline`.
    """
    logger.info(
        "Local execution: %d stage(s), max_workers=%d",
        len(stages),
        config.execution.max_workers,
    )

    return run_pipeline(
        config=config,
        stages=stages,
        state=state,
        force=force,
        dry_run=dry_run,
    )
