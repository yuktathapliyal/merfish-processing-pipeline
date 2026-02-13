"""Configuration loader for the three-layer merFISH pipeline config system.

The three layers are merged in order:

1. **Microscope** -- hardware-specific defaults (``configs/microscopes/<name>.yaml``)
2. **Experiment** -- per-experiment settings (user-provided YAML path)
3. **Execution** -- runtime / scheduler profile (``configs/profiles/<profile>.yaml``)

An optional *overrides* dict (typically from CLI flags) is deep-merged on top
before the final ``PipelineConfig`` is validated.
"""

from __future__ import annotations

import copy
import logging
from pathlib import Path
from typing import Any, Optional

import yaml

from merfish_pipeline.config.defaults import MICROSCOPE_DEFAULTS, VALID_MICROSCOPES
from merfish_pipeline.config.schema import (
    ExperimentConfig,
    ExecutionConfig,
    MicroscopeConfig,
    PipelineConfig,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* into a **copy** of *base*.

    - Dict values are merged recursively.
    - All other types in *override* replace the value in *base*.
    """
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _load_yaml(path: Path) -> dict:
    """Load a YAML file and return its contents as a dict.

    Returns an empty dict when the file is empty.
    """
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return data if isinstance(data, dict) else {}


# ---------------------------------------------------------------------------
# Public: find configs directory
# ---------------------------------------------------------------------------

def find_configs_dir() -> Path:
    """Locate the ``configs/`` directory shipped with the package.

    Search order:
    1. ``<package_root>/../../configs`` (source-tree layout where
       *package_root* is ``src/merfish_pipeline``).
    2. Current working directory ``./configs``.

    Raises ``FileNotFoundError`` if no ``configs/`` directory is found.
    """
    # Relative to *this* file: config/loader.py -> merfish_pipeline -> src -> repo root
    package_dir = Path(__file__).resolve().parent.parent  # merfish_pipeline/
    candidates = [
        package_dir.parent.parent / "configs",  # src/../configs  (repo root)
        Path.cwd() / "configs",
    ]
    for candidate in candidates:
        if candidate.is_dir():
            return candidate

    raise FileNotFoundError(
        "Could not locate the 'configs/' directory. Searched: "
        + ", ".join(str(c) for c in candidates)
    )


# ---------------------------------------------------------------------------
# Layer 1 -- Microscope
# ---------------------------------------------------------------------------

def load_microscope_config(microscope_name: str) -> MicroscopeConfig:
    """Load a microscope configuration by name.

    Looks for ``configs/microscopes/<microscope_name>.yaml``.  Falls back to
    the built-in ``MICROSCOPE_DEFAULTS`` when the YAML file is not present.

    Parameters
    ----------
    microscope_name:
        One of ``"oni"``, ``"nikon"``, or ``"andor"``.

    Returns
    -------
    MicroscopeConfig
    """
    if microscope_name not in VALID_MICROSCOPES:
        raise ValueError(
            f"Unknown microscope {microscope_name!r}. "
            f"Valid options: {VALID_MICROSCOPES}"
        )

    try:
        configs_dir = find_configs_dir()
        yaml_path = configs_dir / "microscopes" / f"{microscope_name}.yaml"
        if yaml_path.is_file():
            data = _load_yaml(yaml_path)
            logger.info("Loaded microscope config from %s", yaml_path)
            return MicroscopeConfig.model_validate(data)
    except FileNotFoundError:
        pass

    # Fallback to built-in defaults
    logger.info(
        "Using built-in defaults for microscope %r (no YAML found).",
        microscope_name,
    )
    if microscope_name not in MICROSCOPE_DEFAULTS:
        raise ValueError(
            f"No built-in defaults for microscope {microscope_name!r}."
        )
    return MicroscopeConfig.model_validate(MICROSCOPE_DEFAULTS[microscope_name])


# ---------------------------------------------------------------------------
# Layer 2 -- Experiment
# ---------------------------------------------------------------------------

def load_experiment_config(path: Path) -> ExperimentConfig:
    """Load an experiment configuration from a YAML file.

    Parameters
    ----------
    path:
        Absolute or relative path to the experiment YAML file.

    Returns
    -------
    ExperimentConfig
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Experiment config not found: {path}")
    data = _load_yaml(path)
    logger.info("Loaded experiment config from %s", path)
    return ExperimentConfig.model_validate(data)


# ---------------------------------------------------------------------------
# Layer 3 -- Execution profile
# ---------------------------------------------------------------------------

def load_execution_config(profile: str) -> ExecutionConfig:
    """Load an execution profile by name.

    Looks for ``configs/profiles/<profile>.yaml``.  Returns a default
    ``ExecutionConfig`` when the file does not exist (so that ``"local"``
    always works without an explicit YAML file).

    Parameters
    ----------
    profile:
        Profile name, e.g. ``"local"`` or ``"slurm"``.

    Returns
    -------
    ExecutionConfig
    """
    try:
        configs_dir = find_configs_dir()
        yaml_path = configs_dir / "profiles" / f"{profile}.yaml"
        if yaml_path.is_file():
            data = _load_yaml(yaml_path)
            logger.info("Loaded execution profile from %s", yaml_path)
            return ExecutionConfig.model_validate(data)
    except FileNotFoundError:
        pass

    if profile == "local":
        logger.info("Using default local execution profile.")
        return ExecutionConfig()

    raise FileNotFoundError(
        f"Execution profile {profile!r} not found and no built-in default exists."
    )


# ---------------------------------------------------------------------------
# Merged -- PipelineConfig
# ---------------------------------------------------------------------------

def load_pipeline_config(
    experiment_path: Path,
    profile: str = "local",
    overrides: Optional[dict[str, Any]] = None,
) -> PipelineConfig:
    """Load and merge all three configuration layers into a ``PipelineConfig``.

    Merge order (later wins):
    1. Microscope defaults (determined by ``experiment.experiment.microscope``)
    2. Experiment YAML
    3. Execution profile YAML
    4. CLI *overrides* dict

    Parameters
    ----------
    experiment_path:
        Path to the experiment YAML file (Layer 2).
    profile:
        Name of the execution profile (Layer 3).  Defaults to ``"local"``.
    overrides:
        Optional dict of CLI / programmatic overrides that are deep-merged on
        top of the combined configuration before validation.

    Returns
    -------
    PipelineConfig
    """
    # --- Layer 2: Experiment ---
    experiment_cfg = load_experiment_config(experiment_path)
    microscope_name = experiment_cfg.experiment.microscope

    # --- Layer 1: Microscope ---
    microscope_cfg = load_microscope_config(microscope_name)

    # --- Layer 3: Execution ---
    execution_cfg = load_execution_config(profile)

    # Build the merged dict.  We start from the experiment data, then inject
    # the resolved microscope and execution objects.
    merged: dict[str, Any] = experiment_cfg.model_dump(mode="python")
    merged["microscope"] = microscope_cfg.model_dump(mode="python")
    merged["execution"] = execution_cfg.model_dump(mode="python")["execution"]

    # --- CLI overrides ---
    if overrides:
        merged = _deep_merge(merged, overrides)

    logger.info(
        "Building PipelineConfig (microscope=%s, profile=%s).",
        microscope_name,
        profile,
    )
    return PipelineConfig.model_validate(merged)
