"""Pydantic v2 models for the three-layer merFISH pipeline configuration system.

Layer 1 -- MicroscopeConfig   (configs/microscopes/*.yaml)
Layer 2 -- ExperimentConfig   (configs/experiments/*.yaml)
Layer 3 -- ExecutionConfig    (configs/profiles/*.yaml)
Merged  -- PipelineConfig     (all three layers combined)
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, field_validator, model_validator

from merfish_pipeline.config.defaults import (
    DEFAULT_DERIVED_PATHS,
    VALID_STAGES,
)


# ============================================================================
# Layer 1 -- Microscope
# ============================================================================

class MicroscopeConfig(BaseModel):
    """Hardware-specific settings loaded from ``configs/microscopes/<name>.yaml``."""

    model_config = {"extra": "forbid"}

    file_pattern: str
    position_file_pattern: Optional[str] = None
    position_format: Literal["csv", "xlsx", "ims"] = "csv"
    stage_x_heading: Optional[str] = "stage_pos_x"
    stage_y_heading: Optional[str] = "stage_pos_y"
    log_file_name: str
    microns_per_pixel: float
    image_dimensions: list[int]
    flip_horizontal: bool = False
    flip_vertical: bool = False
    transpose: bool = False
    default_bead_channel: str = "488nm, Raw"
    default_data_org_template: str
    default_analysis_template: str
    default_microscope_template: str


# ============================================================================
# Layer 2 -- Experiment (nested sub-models first)
# ============================================================================

class AndorConfig(BaseModel):
    """Andor-specific raw-data handling settings."""

    model_config = {"extra": "forbid"}

    channel_order: list[int] = [0, 2, 1]
    round_patterns: list[str] = [
        r"(\d+)(?:st|nd|rd|th)\s*round",
        r"round\s*(\d+)",
        r"R(\d+)",
    ]


class ExperimentInfo(BaseModel):
    model_config = {"extra": "forbid"}

    name: str
    microscope: Literal["oni", "nikon", "andor"]


class PathsConfig(BaseModel):
    model_config = {"extra": "forbid"}

    raw_data_dir: Path
    output_dir: Path
    remapped_data_dir: Optional[Path] = None
    merlin_data_dir: Optional[Path] = None
    parameters_dir: Optional[Path] = None


class RawDataConfig(BaseModel):
    model_config = {"extra": "forbid"}

    bead_channel_folder: str
    data_org_template: Optional[Path] = None
    stage_file: Optional[Path] = None
    stage_x_heading: Optional[str] = None
    stage_y_heading: Optional[str] = None
    andor: Optional[AndorConfig] = None


class MerlinConfig(BaseModel):
    model_config = {"extra": "forbid"}

    codebook_template: Optional[Path] = None
    analysis_template: Optional[Path] = None
    microscope_template: Optional[Path] = None
    cores: int = 100


class FocusQCConfig(BaseModel):
    model_config = {"extra": "forbid"}

    enabled: bool = True
    sigma: float = 1.0
    ksize: int = 3


class StitchConfig(BaseModel):
    model_config = {"extra": "forbid"}

    enabled: bool = True
    group_by: Literal["ir", "z"] = "ir"
    ir_range: Optional[list[int]] = None
    z_range: Optional[list[int]] = None
    images_dir: Optional[Path] = None
    position_file: Optional[Path] = None


class InspectPositionsConfig(BaseModel):
    model_config = {"extra": "forbid"}

    enabled: bool = True
    log_file: Optional[Path] = None
    rounds_to_check: Optional[list[int]] = None
    trajectory_z_slices: Optional[list[int]] = None  # default: first 3 z-slices


class ReregistrationConfig(BaseModel):
    model_config = {"extra": "forbid"}

    enabled: bool = False
    total_z: Optional[int] = None
    target_z: Optional[int] = None
    strict_missing: bool = False
    dry_run: bool = False


class FilterBarcodesConfig(BaseModel):
    model_config = {"extra": "forbid"}

    enabled: bool = False
    mode: Literal["any", "all"] = "any"
    barcodes_file: Optional[Path] = None
    zmap_file: Optional[Path] = None


class SegmentationConfig(BaseModel):
    model_config = {"extra": "forbid"}

    enabled: bool = False
    aligned_images_dir: Optional[Path] = None  # required when enabled
    nuclei_bit: int = 0
    total_bits: int = 16
    exclude_bits: list[int] = []  # additional bits to exclude from cytoplasm sum (0-indexed)
    median_kernel: int = 3
    model_type: str = "cpsam"
    diameter: Optional[int] = None
    batch_size: int = 8
    stitch_threshold: float = 0.5
    # "3d" = segment all z-slices with 2D+stitch (default)
    # "2d" = segment only reference_z_slice, output a single 2D mask
    mode: Literal["2d", "3d"] = "3d"
    reference_z_slice: Optional[int] = None  # required when mode="2d" (see z_indexing)
    z_indexing: Literal[0, 1] = 1  # 0 = 0-indexed, 1 = 1-indexed (default)

    @field_validator("median_kernel")
    @classmethod
    def _check_median_kernel(cls, v: int) -> int:
        if v <= 0 or v % 2 == 0:
            raise ValueError(
                f"median_kernel must be a positive odd integer, got {v}"
            )
        return v


class CorrelationConfig(BaseModel):
    model_config = {"extra": "forbid"}

    enabled: bool = False
    bulk_file: Optional[Path] = None
    distance_thresholds: list[float] = [0.5167, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25]


class PipelineStagesConfig(BaseModel):
    model_config = {"extra": "forbid"}

    stages: list[str] = []
    dry_run: bool = False
    force: bool = False

    @model_validator(mode="after")
    def _validate_stage_names(self) -> "PipelineStagesConfig":
        for name in self.stages:
            if name not in VALID_STAGES:
                raise ValueError(
                    f"Unknown pipeline stage {name!r}. "
                    f"Valid stages: {VALID_STAGES}"
                )
        return self


class ExperimentConfig(BaseModel):
    """Full experiment-level configuration (Layer 2)."""

    model_config = {"extra": "forbid"}

    experiment: ExperimentInfo
    paths: PathsConfig
    raw_data: RawDataConfig
    merlin: MerlinConfig = MerlinConfig()
    focus_qc: FocusQCConfig = FocusQCConfig()
    stitch: StitchConfig = StitchConfig()
    inspect_positions: InspectPositionsConfig = InspectPositionsConfig()
    reregistration: ReregistrationConfig = ReregistrationConfig()
    filter_barcodes: FilterBarcodesConfig = FilterBarcodesConfig()
    segmentation: SegmentationConfig = SegmentationConfig()
    correlation: CorrelationConfig = CorrelationConfig()
    pipeline: PipelineStagesConfig = PipelineStagesConfig()


# ============================================================================
# Layer 3 -- Execution
# ============================================================================

class SlurmConfig(BaseModel):
    """SLURM scheduler settings."""

    model_config = {"extra": "forbid"}

    partition: str
    time: str
    mem_per_cpu: str
    cpus_per_task: int


class ExecutionSettings(BaseModel):
    model_config = {"extra": "forbid"}

    mode: Literal["local", "slurm"] = "local"
    max_workers: int = 8
    slurm: Optional[SlurmConfig] = None
    scratch_dir: Optional[Path] = None


class ExecutionConfig(BaseModel):
    """Execution profile configuration (Layer 3)."""

    model_config = {"extra": "forbid"}

    execution: ExecutionSettings = ExecutionSettings()


# ============================================================================
# Merged -- PipelineConfig
# ============================================================================

class PipelineConfig(BaseModel):
    """Fully resolved pipeline configuration that merges all three layers.

    After validation the following guarantees hold:
    - ``paths.remapped_data_dir``, ``paths.merlin_data_dir``,
      ``paths.parameters_dir`` and ``logs_dir`` are all set (derived from
      ``paths.output_dir`` when not explicitly provided).
    - Microscope defaults for template files are applied when the experiment
      config does not specify them.
    """

    model_config = {"extra": "forbid"}

    # Layer 2 fields
    experiment: ExperimentInfo
    paths: PathsConfig
    raw_data: RawDataConfig
    merlin: MerlinConfig = MerlinConfig()
    focus_qc: FocusQCConfig = FocusQCConfig()
    stitch: StitchConfig = StitchConfig()
    inspect_positions: InspectPositionsConfig = InspectPositionsConfig()
    reregistration: ReregistrationConfig = ReregistrationConfig()
    filter_barcodes: FilterBarcodesConfig = FilterBarcodesConfig()
    segmentation: SegmentationConfig = SegmentationConfig()
    correlation: CorrelationConfig = CorrelationConfig()
    pipeline: PipelineStagesConfig = PipelineStagesConfig()

    # Layer 1 -- resolved microscope settings
    microscope: MicroscopeConfig

    # Layer 3 -- execution settings
    execution: ExecutionSettings = ExecutionSettings()

    # Derived
    logs_dir: Optional[Path] = None

    @model_validator(mode="after")
    def _resolve_derived_paths(self) -> "PipelineConfig":
        """Fill in derived paths from *output_dir* when they are ``None``."""
        output = self.paths.output_dir

        if self.paths.remapped_data_dir is None:
            self.paths.remapped_data_dir = output / DEFAULT_DERIVED_PATHS["remapped_data_dir"]
        if self.paths.merlin_data_dir is None:
            self.paths.merlin_data_dir = output / DEFAULT_DERIVED_PATHS["merlin_data_dir"]
        if self.paths.parameters_dir is None:
            self.paths.parameters_dir = output / DEFAULT_DERIVED_PATHS["parameters_dir"]
        if self.logs_dir is None:
            self.logs_dir = output / DEFAULT_DERIVED_PATHS["logs_dir"]

        return self

    @model_validator(mode="after")
    def _apply_microscope_defaults(self) -> "PipelineConfig":
        """Apply microscope-level template defaults when the experiment does
        not specify them explicitly."""
        mic = self.microscope

        if self.raw_data.data_org_template is None:
            self.raw_data.data_org_template = Path(mic.default_data_org_template)
        if self.merlin.analysis_template is None:
            self.merlin.analysis_template = Path(mic.default_analysis_template)
        if self.merlin.microscope_template is None:
            self.merlin.microscope_template = Path(mic.default_microscope_template)

        # Stage heading defaults from microscope when experiment leaves them unset.
        if self.raw_data.stage_x_heading is None:
            self.raw_data.stage_x_heading = mic.stage_x_heading
        if self.raw_data.stage_y_heading is None:
            self.raw_data.stage_y_heading = mic.stage_y_heading

        return self
