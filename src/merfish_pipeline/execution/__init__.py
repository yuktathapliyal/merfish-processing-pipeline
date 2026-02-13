"""Execution engine package.

Submodules
----------
state
    Pipeline state tracking (JSON persistence).
runner
    Reusable pipeline orchestration logic.
local
    Local (in-process) execution backend.
slurm
    SLURM job generation and submission script builder.
"""

from merfish_pipeline.execution.local import run_local
from merfish_pipeline.execution.runner import resolve_stages, run_pipeline
from merfish_pipeline.execution.slurm import generate_slurm_script
from merfish_pipeline.execution.state import PipelineState

__all__ = [
    "PipelineState",
    "resolve_stages",
    "run_pipeline",
    "run_local",
    "generate_slurm_script",
]
