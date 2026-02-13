"""Stage registry — maps stage names to their implementing classes."""

from __future__ import annotations

from typing import Type

STAGE_REGISTRY: dict[str, Type] = {}


def register_stage(name: str):
    """Class decorator that registers a pipeline stage under *name*.

    Usage::

        @register_stage("focus_qc")
        class FocusQCStage(PipelineStage):
            description = "Per-FOV focus quality check"
            ...

    The decorator sets the class-level ``name`` attribute and adds the class
    to :data:`STAGE_REGISTRY`.
    """

    def decorator(cls):
        if name in STAGE_REGISTRY:
            raise ValueError(
                f"Stage name '{name}' is already registered to "
                f"{STAGE_REGISTRY[name].__name__}"
            )
        STAGE_REGISTRY[name] = cls
        cls.name = name
        return cls

    return decorator


def get_stage(name: str) -> Type:
    """Return the stage class registered under *name*.

    Raises
    ------
    KeyError
        If no stage is registered with the given name.
    """
    try:
        return STAGE_REGISTRY[name]
    except KeyError:
        available = ", ".join(sorted(STAGE_REGISTRY)) or "(none)"
        raise KeyError(
            f"No stage registered with name '{name}'. Available stages: {available}"
        ) from None


def list_stages() -> list[str]:
    """Return a sorted list of all registered stage names."""
    return sorted(STAGE_REGISTRY)
