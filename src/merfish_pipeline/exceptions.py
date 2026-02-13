"""Custom exception hierarchy for the merFISH processing pipeline."""


class MerfishPipelineError(Exception):
    """Base exception for all pipeline errors."""


class ConfigError(MerfishPipelineError):
    """Raised when pipeline configuration is invalid or missing."""


class InputError(MerfishPipelineError):
    """Raised when required input files are missing or malformed."""


class StageError(MerfishPipelineError):
    """Raised when a pipeline stage fails during execution.

    Attributes:
        stage_name: Name of the stage that failed.
        message: Human-readable description of the failure.
        original: The underlying exception, if any.
    """

    def __init__(self, stage_name: str, message: str, original: Exception = None):
        self.stage_name = stage_name
        self.message = message
        self.original = original
        detail = f"Stage '{stage_name}' failed: {message}"
        if original is not None:
            detail += f" (caused by {type(original).__name__}: {original})"
        super().__init__(detail)


class SlurmError(MerfishPipelineError):
    """Raised when a SLURM job submission or execution fails."""
