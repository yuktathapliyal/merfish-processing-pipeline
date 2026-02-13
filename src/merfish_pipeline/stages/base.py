"""Abstract base class and result dataclass for pipeline stages."""

from __future__ import annotations

import hashlib
import json
import logging
import platform
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any


@dataclass
class StageResult:
    """Outcome of a single stage execution.

    Attributes:
        status: One of ``"completed"``, ``"failed"``, or ``"skipped"``.
        output_files: Paths (as strings) to files produced by the stage.
        metadata: Arbitrary key/value pairs recorded by the stage.
        error: Error message if the stage failed, empty string otherwise.
    """

    status: str  # "completed", "failed", "skipped"
    output_files: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    error: str = ""


class PipelineStage(ABC):
    """Abstract base class that every pipeline stage must subclass.

    Subclasses are expected to define ``name`` and ``description`` as class
    attributes (typically set via the ``@register_stage`` decorator for
    ``name``).

    Parameters
    ----------
    config:
        A Pydantic configuration model (must support ``model_dump()``).
    state:
        Optional :class:`~merfish_pipeline.execution.state.PipelineState`
        instance for cross-stage state tracking.
    """

    name: str = ""
    description: str = ""

    def __init__(self, config: Any, state: Any = None):
        self.config = config
        self.state = state
        self.logger = logging.getLogger(f"merfish_pipeline.stages.{self.name}")

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    @abstractmethod
    def validate_inputs(self) -> list[str]:
        """Return a list of error messages.  Empty list means inputs are valid."""

    @abstractmethod
    def check_outputs_exist(self) -> bool:
        """Return ``True`` if all expected outputs already exist (skip logic)."""

    @abstractmethod
    def run(self, dry_run: bool = False) -> StageResult:
        """Execute the stage and return a :class:`StageResult`."""

    # ------------------------------------------------------------------
    # Shared helpers
    # ------------------------------------------------------------------

    def should_skip(self, force: bool = False) -> bool:
        """Decide whether this stage can be skipped.

        Returns ``False`` when *force* is ``True``; otherwise delegates to
        :meth:`check_outputs_exist`.
        """
        if force:
            return False
        return self.check_outputs_exist()

    def get_output_dir(self) -> Path:
        """Return the output directory for this stage."""
        return Path(self.config.paths.output_dir) / self.name

    def write_run_metadata(
        self,
        result: StageResult,
        start_time: datetime,
        parameters: dict[str, Any] | None = None,
    ) -> None:
        """Persist run metadata as JSON inside the stage output directory.

        The file is written to ``<output_dir>/<stage_name>/run_metadata.json``
        and contains timing information, status, config hash, output file
        list, and the host on which the stage ran.
        """
        end_time = datetime.now()
        runtime_seconds = round((end_time - start_time).total_seconds(), 3)

        config_dict = self.config.model_dump()
        config_json = json.dumps(config_dict, sort_keys=True, default=str)
        config_hash = hashlib.sha256(config_json.encode("utf-8")).hexdigest()

        metadata = {
            "stage_name": self.name,
            "started_at": start_time.isoformat(),
            "completed_at": end_time.isoformat(),
            "status": result.status,
            "pipeline_version": _get_pipeline_version(),
            "config_hash": config_hash,
            "output_files": result.output_files,
            "parameters": parameters or {},
            "runtime_seconds": runtime_seconds,
            "hostname": platform.node(),
        }

        output_dir = self.get_output_dir()
        output_dir.mkdir(parents=True, exist_ok=True)

        metadata_path = output_dir / "run_metadata.json"
        metadata_path.write_text(
            json.dumps(metadata, indent=2, default=str) + "\n",
            encoding="utf-8",
        )
        self.logger.debug("Wrote run metadata to %s", metadata_path)


def _get_pipeline_version() -> str:
    """Return the installed pipeline version, falling back to 'unknown'."""
    try:
        from merfish_pipeline import __version__

        return __version__
    except Exception:
        return "unknown"
