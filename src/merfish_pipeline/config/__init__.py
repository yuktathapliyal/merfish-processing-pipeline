"""Configuration system for the merFISH processing pipeline.

Three-layer config:
  Layer 1 -- MicroscopeConfig   (hardware defaults)
  Layer 2 -- ExperimentConfig   (per-experiment settings)
  Layer 3 -- ExecutionConfig    (runtime / scheduler profile)
  Merged  -- PipelineConfig     (fully resolved)
"""

from merfish_pipeline.config.schema import (
    AndorConfig,
    CorrelationConfig,
    ExecutionConfig,
    ExecutionSettings,
    ExperimentConfig,
    ExperimentInfo,
    FilterBarcodesConfig,
    FocusQCConfig,
    InspectPositionsConfig,
    MerlinConfig,
    MicroscopeConfig,
    PathsConfig,
    PipelineConfig,
    PipelineStagesConfig,
    RawDataConfig,
    ReregistrationConfig,
    SegmentationConfig,
    SlurmConfig,
    StitchConfig,
)
from merfish_pipeline.config.loader import (
    find_configs_dir,
    load_experiment_config,
    load_execution_config,
    load_microscope_config,
    load_pipeline_config,
)

__all__ = [
    # Schema models
    "AndorConfig",
    "CorrelationConfig",
    "ExecutionConfig",
    "ExecutionSettings",
    "ExperimentConfig",
    "ExperimentInfo",
    "FilterBarcodesConfig",
    "FocusQCConfig",
    "InspectPositionsConfig",
    "MerlinConfig",
    "MicroscopeConfig",
    "PathsConfig",
    "PipelineConfig",
    "PipelineStagesConfig",
    "RawDataConfig",
    "ReregistrationConfig",
    "SegmentationConfig",
    "SlurmConfig",
    "StitchConfig",
    # Loader functions
    "find_configs_dir",
    "load_experiment_config",
    "load_execution_config",
    "load_microscope_config",
    "load_pipeline_config",
]
