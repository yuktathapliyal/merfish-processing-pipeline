"""Pipeline state tracking — persists stage execution status to a JSON file."""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

from merfish_pipeline.stages.base import StageResult

logger = logging.getLogger("merfish_pipeline.execution.state")

STATE_FILENAME = ".pipeline_state.json"


class PipelineState:
    """Track per-stage execution status in a JSON file on disk.

    Parameters
    ----------
    state_path:
        Absolute path to the JSON state file.
    """

    def __init__(self, state_path: Path):
        self.state_path = Path(state_path)
        self._data: dict[str, Any] = {"experiment": "", "stages": {}}

        if self.state_path.exists():
            self._load()

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def load_or_create(cls, output_dir: Path) -> PipelineState:
        """Load an existing state file from *output_dir* or create a new one.

        The state file is stored as ``{output_dir}/.pipeline_state.json``.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        state_path = output_dir / STATE_FILENAME
        return cls(state_path)

    # ------------------------------------------------------------------
    # Stage lifecycle
    # ------------------------------------------------------------------

    def mark_started(self, stage_name: str) -> None:
        """Record that *stage_name* has begun execution."""
        self._data["stages"][stage_name] = {
            "status": "running",
            "started_at": datetime.now().isoformat(),
            "completed_at": "",
            "output_files": [],
            "metadata": {},
            "error": "",
        }
        self.save()

    def mark_completed(self, stage_name: str, result: StageResult) -> None:
        """Record that *stage_name* finished successfully."""
        stage = self._data["stages"].setdefault(stage_name, {})
        stage["status"] = "completed"
        stage["completed_at"] = datetime.now().isoformat()
        stage["output_files"] = result.output_files
        stage["metadata"] = result.metadata
        stage["error"] = ""
        self.save()

    def mark_failed(self, stage_name: str, error: str) -> None:
        """Record that *stage_name* failed with *error*."""
        stage = self._data["stages"].setdefault(stage_name, {})
        stage["status"] = "failed"
        stage["completed_at"] = datetime.now().isoformat()
        stage["error"] = error
        self.save()

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def is_completed(self, stage_name: str) -> bool:
        """Return ``True`` if *stage_name* has status ``"completed"``."""
        stage = self._data["stages"].get(stage_name, {})
        return stage.get("status") == "completed"

    def get_stage_status(self, stage_name: str) -> dict[str, Any]:
        """Return the full status dict for *stage_name*, or empty dict."""
        return self._data["stages"].get(stage_name, {})

    def summary(self) -> dict[str, Any]:
        """Return a summary of pipeline execution state.

        Includes the experiment name, per-stage status strings, and aggregate
        counts of completed / failed / running stages.
        """
        stages = self._data.get("stages", {})
        status_counts: dict[str, int] = {}
        stage_statuses: dict[str, str] = {}
        for name, info in stages.items():
            status = info.get("status", "unknown")
            stage_statuses[name] = status
            status_counts[status] = status_counts.get(status, 0) + 1

        return {
            "experiment": self._data.get("experiment", ""),
            "stages": stage_statuses,
            "counts": status_counts,
        }

    # ------------------------------------------------------------------
    # Property access
    # ------------------------------------------------------------------

    @property
    def experiment(self) -> str:
        return self._data.get("experiment", "")

    @experiment.setter
    def experiment(self, value: str) -> None:
        self._data["experiment"] = value

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self) -> None:
        """Write the current state to disk."""
        self.state_path.parent.mkdir(parents=True, exist_ok=True)
        self.state_path.write_text(
            json.dumps(self._data, indent=2, default=str) + "\n",
            encoding="utf-8",
        )
        logger.debug("State saved to %s", self.state_path)

    def _load(self) -> None:
        """Read state from disk."""
        try:
            text = self.state_path.read_text(encoding="utf-8")
            self._data = json.loads(text)
        except (json.JSONDecodeError, OSError) as exc:
            logger.warning("Could not load state from %s: %s", self.state_path, exc)
            self._data = {"experiment": "", "stages": {}}
