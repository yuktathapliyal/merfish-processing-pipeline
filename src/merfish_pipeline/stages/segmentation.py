"""``segmentation`` stage -- Cellpose cell segmentation on microscopy images.

This stage combines preprocessing (extracting nuclei/cytoplasm channels from
aligned MERlin TIFF stacks) and 2D+stitch Cellpose segmentation into a single
pipeline stage.

Algorithm
---------
**Preprocessing** (per FOV):

1. Load an aligned multi-frame TIFF stack produced by MERlin.
2. Extract the nuclei channel from the specified bit/frame index.
3. Build a cytoplasm channel by summing all remaining bits across Z.
4. Apply a median filter (configurable kernel size) to both channels.
5. Normalize each channel to [0, 1] via min-max normalization.
6. Output a 4-D volume ``(Z, C=2, Y, X)`` where ``C[0]`` = cytoplasm and
   ``C[1]`` = nuclei.

**Segmentation** (per FOV):

1. Initialize a Cellpose model (model type from config, GPU preferred).
2. Run Cellpose in 2D+stitch mode: segment each z-slice independently,
   then stitch masks across z using an overlap threshold.
3. Save the resulting ``uint16`` mask TIFF alongside the preprocessed volume.

**Single-slice mode** (optional):

When ``reference_z_slice`` is set in config, only that z-slice is preprocessed
and segmented.  The output mask is 2-D ``(Y, X)`` instead of 3-D ``(Z, Y, X)``.
This is useful when only one z-slice has sufficient image quality (e.g. Nikon).
Downstream cell assignment applies the single mask to barcodes from all z-slices.

Outputs
-------
- ``{output_dir}/segmentation/preprocessed/``  -- preprocessed volumes.
- ``{output_dir}/segmentation/masks/``         -- per-FOV segmentation masks
  (3-D if all z-slices, 2-D if single-slice mode).
- ``{output_dir}/segmentation/run_metadata.json`` -- timing and parameters.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

try:
    import cv2
except ImportError:
    cv2 = None  # checked in validate_inputs()

from merfish_pipeline.io.tiff_io import read_tiff, write_tiff
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helper: find aligned images
# ---------------------------------------------------------------------------


def _find_aligned_images(input_dir: Path, pattern: str) -> list[Path]:
    """Discover aligned-image TIFF stacks in *input_dir*.

    Parameters
    ----------
    input_dir:
        Directory to search (typically a MERlin output directory).
    pattern:
        Glob pattern for the aligned image files, e.g.
        ``"aligned_images*.tif"`` or ``"aligned_images*.tiff"``.

    Returns
    -------
    list[Path]
        Sorted list of matching paths.
    """
    if not input_dir.is_dir():
        return []
    matches = sorted(input_dir.glob(pattern), key=lambda p: p.name)
    return matches


# ---------------------------------------------------------------------------
# Helper: normalization
# ---------------------------------------------------------------------------


def _normalize(image: np.ndarray) -> np.ndarray:
    """Min-max normalize *image* to the range [0, 1] as float32.

    Parameters
    ----------
    image:
        Input array of any numeric dtype.

    Returns
    -------
    np.ndarray
        Float32 array with values in [0, 1].
    """
    img = image.astype(np.float32)
    lo = img.min()
    hi = img.max()
    return (img - lo) / (hi - lo + 1e-8)


# ---------------------------------------------------------------------------
# Helper: preprocess a single volume
# ---------------------------------------------------------------------------


def _preprocess_volume(
    tiff_path: Path,
    nuclei_bit: int,
    total_bits: int,
    mf_kernel: int,
    z_index: int | None = None,
    exclude_bits: list[int] | None = None,
) -> np.ndarray:
    """Extract and preprocess nuclei + cytoplasm channels from one FOV stack.

    Matches the reference ``processZstacksForCellpose.py`` algorithm:
    each plane is median-filtered and normalized independently (per-slice,
    per-bit), and the cytoplasm channel is the sum of individually
    normalized non-nuclei bit-planes.

    Parameters
    ----------
    tiff_path:
        Path to the aligned-image TIFF stack.
    nuclei_bit:
        0-based index of the nuclei channel within the bit dimension.
    total_bits:
        Expected total number of bits (channels) in the stack.
    mf_kernel:
        Kernel size for :func:`cv2.medianBlur`.  Must be a positive odd
        integer.
    z_index:
        If provided (0-based), extract only this z-slice.  The output will
        still be 4-D with ``Z=1``.
    exclude_bits:
        Additional 0-based bit indices to exclude from the cytoplasm sum
        (e.g. fiducial beads).  ``nuclei_bit`` is always excluded.

    Returns
    -------
    np.ndarray
        4-D float32 array with shape ``(Z, 2, Y, X)``.
        Channel 0 is the preprocessed cytoplasm signal and channel 1 is
        the preprocessed nuclei signal.  When *z_index* is set, ``Z=1``.

    Raises
    ------
    ValueError
        If the TIFF shape is incompatible with the expected layout.
    """
    if cv2 is None:
        raise RuntimeError(
            "OpenCV (cv2) is required for segmentation preprocessing. "
            "Install it with: pip install opencv-python-headless"
        )

    if exclude_bits is None:
        exclude_bits = []

    raw = read_tiff(tiff_path)
    logger.debug("Loaded %s with shape %s dtype %s", tiff_path.name, raw.shape, raw.dtype)

    # Determine Z, bits, Y, X from the raw shape.
    # MERlin writes aligned_images as bit-major, z-minor:
    #   bit0_z0, bit0_z1, ..., bit0_zN, bit1_z0, ..., bitM_zN
    # So the flat 3-D shape is (bits*Z, Y, X) with bits changing slowest.
    # We reshape to (bits, Z, Y, X) then transpose to (Z, bits, Y, X).
    if raw.ndim == 4:
        # MERlin writes 4-D as (bits, Z, Y, X); transpose to (Z, bits, Y, X)
        n_bits, n_z, h, w = raw.shape
        raw = raw.transpose(1, 0, 2, 3)
    elif raw.ndim == 3:
        n_frames, h, w = raw.shape
        if total_bits <= 0:
            raise ValueError(
                f"Cannot determine Z dimension: total_bits must be > 0 (got {total_bits})"
            )
        if n_frames % total_bits != 0:
            raise ValueError(
                f"Frame count {n_frames} is not divisible by total_bits {total_bits} "
                f"for {tiff_path.name}"
            )
        n_z = n_frames // total_bits
        n_bits = total_bits
        raw = raw.reshape(n_bits, n_z, h, w).transpose(1, 0, 2, 3)
    else:
        raise ValueError(
            f"Unexpected TIFF shape {raw.shape} for {tiff_path.name}; "
            f"expected 3-D or 4-D array."
        )

    if nuclei_bit < 0 or nuclei_bit >= n_bits:
        raise ValueError(
            f"nuclei_bit={nuclei_bit} is out of range for stack with {n_bits} bits"
        )

    # Optionally select a single z-slice
    if z_index is not None:
        if z_index < 0 or z_index >= n_z:
            raise ValueError(
                f"reference_z_slice maps to z_index={z_index} which is out of "
                f"range [0, {n_z}) for {tiff_path.name}"
            )
        raw = raw[z_index : z_index + 1, :, :, :]  # keep 4-D: (1, bits, Y, X)
        n_z = 1

    # Determine which bits contribute to the cytoplasm channel
    skip = {nuclei_bit} | set(exclude_bits)
    cyto_bits = [b for b in range(n_bits) if b not in skip]

    # Build output volume with per-slice, per-bit normalization
    # (matches reference processZstacksForCellpose.py)
    vol = np.zeros((n_z, 2, h, w), dtype=np.float32)

    for z in range(n_z):
        # Nuclei channel (channel 1): median filter + normalize this slice
        nuc_plane = raw[z, nuclei_bit].astype(np.float32)
        nuc_plane = cv2.medianBlur(nuc_plane, mf_kernel)
        vol[z, 1] = _normalize(nuc_plane)

        # Cytoplasm channel (channel 0): sum of per-plane normalized bits
        accum = np.zeros((h, w), dtype=np.float32)
        for b in cyto_bits:
            plane = raw[z, b].astype(np.float32)
            plane = cv2.medianBlur(plane, mf_kernel)
            accum += _normalize(plane)
        vol[z, 0] = accum

    return vol


# ---------------------------------------------------------------------------
# Helper: segment a single volume
# ---------------------------------------------------------------------------


def _segment_volume(
    volume: np.ndarray,
    model: Any,
    diameter: int | None,
    batch_size: int,
    stitch_threshold: float,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
) -> np.ndarray:
    """Run Cellpose 2D+stitch segmentation on a preprocessed volume.

    Parameters
    ----------
    volume:
        4-D float32 array ``(Z, C=2, Y, X)`` from :func:`_preprocess_volume`.
    model:
        An initialized ``cellpose.models.CellposeModel`` instance.
    diameter:
        Expected cell diameter in pixels.  ``None`` triggers Cellpose
        auto-estimation.
    batch_size:
        Number of images processed per GPU batch.
    stitch_threshold:
        IoU overlap threshold used to stitch masks across z-slices.
    flow_threshold:
        Cellpose flow error threshold.  Lower values produce fewer,
        higher-confidence cells.
    cellprob_threshold:
        Minimum cell probability to accept a pixel as belonging to a cell.

    Returns
    -------
    np.ndarray
        3-D ``uint16`` mask array ``(Z, Y, X)`` where each unique non-zero
        value represents one segmented cell.
    """
    # Cellpose v4+ requires explicit axis parameters for 4-D input.
    # Our volume shape is (Z, C=2, Y, X): Z=axis 0, C=axis 1.
    masks, flows, styles = model.eval(
        volume,
        z_axis=0,
        channel_axis=1,
        diameter=diameter,
        do_3D=False,
        stitch_threshold=stitch_threshold,
        batch_size=batch_size,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )

    # Ensure output is uint16 for compact storage.
    masks = np.asarray(masks, dtype=np.uint16)
    return masks


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("segmentation")
class SegmentationStage(PipelineStage):
    """Cellpose cell segmentation on aligned MERlin microscopy images."""

    description = "Run Cellpose cell segmentation on aligned microscopy images"

    _DEFAULT_ALIGNED_PATTERN: str = "aligned_images*.tif*"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that aligned-image inputs exist and are accessible."""
        errors: list[str] = []

        seg_cfg = self.config.segmentation

        # --- Config consistency checks ---
        if seg_cfg.mode == "2d" and seg_cfg.reference_z_slice is None:
            errors.append(
                "segmentation.mode is '2d' but reference_z_slice is not set. "
                "Specify which z-slice to segment."
            )

        if seg_cfg.nuclei_bit >= seg_cfg.total_bits:
            errors.append(
                f"segmentation.nuclei_bit ({seg_cfg.nuclei_bit}) must be less "
                f"than total_bits ({seg_cfg.total_bits})."
            )

        # --- Check cv2 availability ---
        try:
            import cv2  # noqa: F401
        except ImportError:
            errors.append(
                "OpenCV (cv2) is not installed. "
                "Install it with: pip install opencv-python-headless"
            )

        # --- Check cellpose availability ---
        try:
            import cellpose  # noqa: F401
        except ImportError:
            errors.append(
                "cellpose is not installed. Install it with: "
                "pip install 'merfish-pipeline[segmentation]' or "
                "pip install cellpose"
            )

        # --- Aligned images directory ---
        input_dir = self._resolve_input_dir()
        if input_dir is None:
            auto_path = (
                Path(self.config.paths.output_dir)
                / "merlin_analysis"
                / self.config.experiment.name
                / "FiducialCorrelationWarp"
                / "images"
            )
            errors.append(
                "Could not find aligned images directory. Either:\n"
                f"  - Set segmentation.aligned_images_dir in your config, or\n"
                f"  - Run MERlin first so that {auto_path} exists."
            )
            return errors

        if not input_dir.exists():
            errors.append(f"Aligned images directory does not exist: {input_dir}")
        elif not input_dir.is_dir():
            errors.append(f"Aligned images path is not a directory: {input_dir}")
        else:
            images = _find_aligned_images(input_dir, self._DEFAULT_ALIGNED_PATTERN)
            if not images:
                errors.append(
                    f"No aligned image TIFFs found in {input_dir} "
                    f"matching pattern '{self._DEFAULT_ALIGNED_PATTERN}'"
                )

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if the masks directory already contains at least one TIFF."""
        masks_dir = self.get_output_dir() / "masks"
        if not masks_dir.exists():
            return False
        mask_files = list(masks_dir.glob("*.tif*"))
        return len(mask_files) > 0

    def run(self, dry_run: bool = False) -> StageResult:
        """Execute preprocessing and Cellpose segmentation."""
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        masks_dir = output_dir / "masks"
        preprocessed_dir = output_dir / "preprocessed"

        seg_cfg = self.config.segmentation
        nuclei_bit = seg_cfg.nuclei_bit
        total_bits = seg_cfg.total_bits
        median_kernel = seg_cfg.median_kernel
        model_type = seg_cfg.model_type
        diameter = seg_cfg.diameter
        batch_size = seg_cfg.batch_size
        stitch_threshold = seg_cfg.stitch_threshold
        flow_threshold = seg_cfg.flow_threshold
        cellprob_threshold = seg_cfg.cellprob_threshold

        # Segmentation mode: "2d" = single reference slice, "3d" = all slices
        single_slice_mode = seg_cfg.mode == "2d"
        z_index: int | None = None
        if single_slice_mode:
            if seg_cfg.reference_z_slice is None:
                return StageResult(
                    status="failed",
                    error="mode is '2d' but reference_z_slice is not set. "
                    "Specify which z-slice to segment.",
                )
            z_index = seg_cfg.reference_z_slice - seg_cfg.z_indexing
            self.logger.info(
                "2D mode: segmenting z-slice %d (0-based) only",
                z_index,
            )
        else:
            self.logger.info("3D mode: segmenting all z-slices with 2D+stitch")

        # ----------------------------------------------------------
        # 0. Check optional cellpose dependency
        # ----------------------------------------------------------
        try:
            from cellpose import models  # noqa: F401
        except ImportError:
            msg = (
                "cellpose is not installed. Install it with: "
                "pip install 'merfish-pipeline[segmentation]' or "
                "pip install cellpose"
            )
            self.logger.error(msg)
            return StageResult(status="failed", error=msg)

        # ----------------------------------------------------------
        # 1. Discover aligned images
        # ----------------------------------------------------------
        input_dir = self._resolve_input_dir()
        if input_dir is None:
            return StageResult(
                status="failed",
                error="Could not find aligned images directory. "
                "Set segmentation.aligned_images_dir or run MERlin first.",
            )
        self.logger.info("Scanning for aligned images in %s", input_dir)

        aligned_images = _find_aligned_images(input_dir, self._DEFAULT_ALIGNED_PATTERN)
        if not aligned_images:
            return StageResult(
                status="failed",
                error=f"No aligned image TIFFs found in {input_dir}",
            )

        self.logger.info("Found %d aligned image stack(s).", len(aligned_images))

        if dry_run:
            mode_desc = f"z_slice={z_index} (0-based)" if single_slice_mode else "all z-slices (2D+stitch)"
            self.logger.info(
                "[DRY RUN] Would preprocess and segment %d FOV(s) "
                "using model_type=%s, diameter=%s, nuclei_bit=%d, total_bits=%d, mode=%s",
                len(aligned_images), model_type, diameter, nuclei_bit, total_bits, mode_desc,
            )
            return StageResult(
                status="skipped",
                metadata={
                    "dry_run": True,
                    "n_fovs": len(aligned_images),
                    "nuclei_bit": nuclei_bit,
                    "total_bits": total_bits,
                    "model_type": model_type,
                    "diameter": diameter,
                    "stitch_threshold": stitch_threshold,
                    "flow_threshold": flow_threshold,
                    "cellprob_threshold": cellprob_threshold,
                    "single_slice_mode": single_slice_mode,
                    "reference_z_index": z_index,
                },
            )

        # ----------------------------------------------------------
        # 2. Initialize Cellpose model
        # ----------------------------------------------------------
        self.logger.info(
            "Initializing Cellpose model (type=%s) ...", model_type,
        )
        model = self._init_cellpose_model(models, model_type)

        # ----------------------------------------------------------
        # 3. Create output directories
        # ----------------------------------------------------------
        masks_dir.mkdir(parents=True, exist_ok=True)
        preprocessed_dir.mkdir(parents=True, exist_ok=True)

        # ----------------------------------------------------------
        # 4. Process each FOV
        # ----------------------------------------------------------
        output_files: list[str] = []
        n_processed = 0
        n_failed = 0
        fov_errors: list[str] = []

        for tiff_path in aligned_images:
            fov_name = tiff_path.stem
            self.logger.info("Processing FOV: %s", fov_name)

            # --- 4a. Preprocess ---
            try:
                volume = _preprocess_volume(
                    tiff_path=tiff_path,
                    nuclei_bit=nuclei_bit,
                    total_bits=total_bits,
                    mf_kernel=median_kernel,
                    z_index=z_index,
                    exclude_bits=list(seg_cfg.exclude_bits),
                )
            except Exception as exc:
                msg = f"Preprocessing failed for {fov_name}: {exc}"
                self.logger.warning(msg)
                fov_errors.append(msg)
                n_failed += 1
                continue

            # Save preprocessed volume
            preproc_path = preprocessed_dir / f"{fov_name}_preprocessed.tif"
            write_tiff(volume, preproc_path)
            output_files.append(str(preproc_path))
            self.logger.debug("Saved preprocessed volume: %s", preproc_path)

            # --- 4b. Segment ---
            try:
                masks = _segment_volume(
                    volume=volume,
                    model=model,
                    diameter=diameter,
                    batch_size=batch_size,
                    stitch_threshold=stitch_threshold,
                    flow_threshold=flow_threshold,
                    cellprob_threshold=cellprob_threshold,
                )
            except Exception as exc:
                msg = f"Segmentation failed for {fov_name}: {exc}"
                self.logger.warning(msg)
                fov_errors.append(msg)
                n_failed += 1
                continue

            # In single-slice mode, squeeze to 2D mask (Y, X)
            if single_slice_mode and masks.ndim == 3 and masks.shape[0] == 1:
                masks = masks[0]

            # Save mask
            mask_path = masks_dir / f"{fov_name}_masks.tif"
            write_tiff(masks, mask_path)
            output_files.append(str(mask_path))
            self.logger.info(
                "Saved mask for %s: %d cells detected (shape %s)",
                fov_name,
                int(masks.max()),
                masks.shape,
            )

            n_processed += 1
            if n_processed % 10 == 0:
                self.logger.info(
                    "  Progress: %d / %d FOVs processed ...",
                    n_processed,
                    len(aligned_images),
                )

        # ----------------------------------------------------------
        # 5. Build result
        # ----------------------------------------------------------
        self.logger.info(
            "Segmentation complete: %d succeeded, %d failed out of %d FOVs.",
            n_processed, n_failed, len(aligned_images),
        )

        if n_processed == 0:
            result = StageResult(
                status="failed",
                output_files=output_files,
                error=(
                    f"All {len(aligned_images)} FOVs failed during segmentation. "
                    f"Errors: {'; '.join(fov_errors)}"
                ),
            )
        else:
            status = "completed"
            result = StageResult(
                status=status,
                output_files=output_files,
                metadata={
                    "n_fovs_total": len(aligned_images),
                    "n_fovs_processed": n_processed,
                    "n_fovs_failed": n_failed,
                    "model_type": model_type,
                    "diameter": diameter,
                    "batch_size": batch_size,
                    "stitch_threshold": stitch_threshold,
                    "flow_threshold": flow_threshold,
                    "cellprob_threshold": cellprob_threshold,
                    "single_slice_mode": single_slice_mode,
                    "reference_z_index": z_index,
                    "fov_errors": fov_errors,
                },
            )

        # ----------------------------------------------------------
        # 6. Write run metadata
        # ----------------------------------------------------------
        self.write_run_metadata(
            result,
            start_time,
            parameters={
                "nuclei_bit": nuclei_bit,
                "total_bits": total_bits,
                "median_kernel": median_kernel,
                "model_type": model_type,
                "diameter": diameter,
                "batch_size": batch_size,
                "stitch_threshold": stitch_threshold,
                "flow_threshold": flow_threshold,
                "cellprob_threshold": cellprob_threshold,
                "single_slice_mode": single_slice_mode,
                "reference_z_index": z_index,
                "aligned_pattern": self._DEFAULT_ALIGNED_PATTERN,
            },
        )

        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_input_dir(self) -> Path | None:
        """Resolve the aligned-images directory.

        Priority:
        1. Explicit ``segmentation.aligned_images_dir`` config override.
        2. Auto-detect at ``{output_dir}/merlin_analysis/{xp}/FiducialCorrelationWarp/images``.
        """
        if self.config.segmentation.aligned_images_dir is not None:
            return Path(self.config.segmentation.aligned_images_dir)

        candidate = (
            Path(self.config.paths.output_dir)
            / "merlin_analysis"
            / self.config.experiment.name
            / "FiducialCorrelationWarp"
            / "images"
        )
        if candidate.is_dir():
            self.logger.info("Auto-detected aligned images at %s", candidate)
            return candidate

        return None

    @staticmethod
    def _init_cellpose_model(models_module: Any, model_type: str) -> Any:
        """Initialize a Cellpose model, falling back to CPU if GPU fails.

        Parameters
        ----------
        models_module:
            The ``cellpose.models`` module (already imported by caller).
        model_type:
            Cellpose model type string, e.g. ``"cyto2"`` or ``"cpsam"``.

        Returns
        -------
        A ``CellposeModel`` instance.
        """
        try:
            model = models_module.CellposeModel(model_type=model_type, gpu=True)
            logger.info("Cellpose model initialized with GPU support.")
            return model
        except Exception as exc:
            logger.warning(
                "GPU initialization failed (%s); falling back to CPU.", exc,
            )
            model = models_module.CellposeModel(model_type=model_type, gpu=False)
            logger.info("Cellpose model initialized on CPU.")
            return model
