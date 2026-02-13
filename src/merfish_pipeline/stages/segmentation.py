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

Outputs
-------
- ``{output_dir}/segmentation/preprocessed/``  -- preprocessed 4-D volumes.
- ``{output_dir}/segmentation/masks/``         -- per-FOV segmentation masks.
- ``{output_dir}/segmentation/run_metadata.json`` -- timing and parameters.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

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
) -> np.ndarray:
    """Extract and preprocess nuclei + cytoplasm channels from one FOV stack.

    Parameters
    ----------
    tiff_path:
        Path to the aligned-image TIFF stack.
    nuclei_bit:
        Index of the nuclei channel within the bit dimension.
    total_bits:
        Expected total number of bits (channels) in the stack.  Used to
        validate the input and to determine which bits form the cytoplasm
        channel.
    mf_kernel:
        Kernel size for the :func:`cv2.medianBlur` filter.  Must be a
        positive odd integer.

    Returns
    -------
    np.ndarray
        4-D float32 array with shape ``(Z, 2, Y, X)``.
        Channel 0 is the preprocessed cytoplasm signal and channel 1 is
        the preprocessed nuclei signal.

    Raises
    ------
    ValueError
        If the TIFF shape is incompatible with the expected layout.
    """
    import cv2

    raw = read_tiff(tiff_path)
    logger.debug("Loaded %s with shape %s dtype %s", tiff_path.name, raw.shape, raw.dtype)

    # Determine Z, bits, Y, X from the raw shape.
    # Typical shapes:
    #   (Z, bits, Y, X)  -- 4-D
    #   (Z*bits, Y, X)   -- 3-D, needs reshaping
    if raw.ndim == 4:
        n_z, n_bits, h, w = raw.shape
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
        raw = raw.reshape(n_z, n_bits, h, w)
    else:
        raise ValueError(
            f"Unexpected TIFF shape {raw.shape} for {tiff_path.name}; "
            f"expected 3-D or 4-D array."
        )

    if nuclei_bit < 0 or nuclei_bit >= n_bits:
        raise ValueError(
            f"nuclei_bit={nuclei_bit} is out of range for stack with {n_bits} bits"
        )

    # Extract nuclei channel: (Z, Y, X)
    nuclei = raw[:, nuclei_bit, :, :].astype(np.float32)

    # Build cytoplasm channel: sum of all bits except nuclei across Z -> (Z, Y, X)
    cyto_bits = [b for b in range(n_bits) if b != nuclei_bit]
    if cyto_bits:
        cyto = raw[:, cyto_bits, :, :].astype(np.float32).sum(axis=1)
    else:
        # Degenerate case: only one bit available.
        cyto = np.zeros_like(nuclei)

    # Apply median filter per z-slice
    for z_idx in range(n_z):
        nuclei[z_idx] = cv2.medianBlur(nuclei[z_idx], mf_kernel)
        cyto[z_idx] = cv2.medianBlur(cyto[z_idx], mf_kernel)

    # Normalize
    nuclei = _normalize(nuclei)
    cyto = _normalize(cyto)

    # Assemble output: (Z, C=2, Y, X)
    volume = np.stack([cyto, nuclei], axis=1)
    return volume


# ---------------------------------------------------------------------------
# Helper: segment a single volume
# ---------------------------------------------------------------------------


def _segment_volume(
    volume: np.ndarray,
    model: Any,
    diameter: int | None,
    batch_size: int,
    stitch_threshold: float,
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

    Returns
    -------
    np.ndarray
        3-D ``uint16`` mask array ``(Z, Y, X)`` where each unique non-zero
        value represents one segmented cell.
    """
    # Cellpose eval expects (Z, C, Y, X) for multi-channel, channels=[1, 2]
    # means: channel 1 = cytoplasm (green), channel 2 = nuclei (red).
    masks, flows, styles = model.eval(
        volume,
        channels=[1, 2],
        diameter=diameter,
        do_3D=False,
        stitch_threshold=stitch_threshold,
        batch_size=batch_size,
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

    # Default parameters when not overridden by upstream config.
    _DEFAULT_NUCLEI_BIT: int = 0
    _DEFAULT_TOTAL_BITS: int = 16
    _DEFAULT_MEDIAN_KERNEL: int = 3
    _DEFAULT_ALIGNED_PATTERN: str = "aligned_images*.tif*"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that aligned-image inputs exist and are accessible."""
        errors: list[str] = []

        input_dir = self._resolve_input_dir()
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
        model_type = seg_cfg.model_type
        diameter = seg_cfg.diameter
        batch_size = seg_cfg.batch_size
        stitch_threshold = seg_cfg.stitch_threshold

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
        self.logger.info("Scanning for aligned images in %s", input_dir)

        aligned_images = _find_aligned_images(input_dir, self._DEFAULT_ALIGNED_PATTERN)
        if not aligned_images:
            return StageResult(
                status="failed",
                error=f"No aligned image TIFFs found in {input_dir}",
            )

        self.logger.info("Found %d aligned image stack(s).", len(aligned_images))

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would preprocess and segment %d FOV(s) "
                "using model_type=%s, diameter=%s, stitch_threshold=%s",
                len(aligned_images), model_type, diameter, stitch_threshold,
            )
            return StageResult(
                status="skipped",
                metadata={
                    "dry_run": True,
                    "n_fovs": len(aligned_images),
                    "model_type": model_type,
                    "diameter": diameter,
                    "stitch_threshold": stitch_threshold,
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
                    nuclei_bit=self._DEFAULT_NUCLEI_BIT,
                    total_bits=self._DEFAULT_TOTAL_BITS,
                    mf_kernel=self._DEFAULT_MEDIAN_KERNEL,
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
                )
            except Exception as exc:
                msg = f"Segmentation failed for {fov_name}: {exc}"
                self.logger.warning(msg)
                fov_errors.append(msg)
                n_failed += 1
                continue

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
                "model_type": model_type,
                "diameter": diameter,
                "batch_size": batch_size,
                "stitch_threshold": stitch_threshold,
                "nuclei_bit": self._DEFAULT_NUCLEI_BIT,
                "total_bits": self._DEFAULT_TOTAL_BITS,
                "median_kernel": self._DEFAULT_MEDIAN_KERNEL,
                "aligned_pattern": self._DEFAULT_ALIGNED_PATTERN,
            },
        )

        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_input_dir(self) -> Path:
        """Resolve the directory containing aligned-image TIFF stacks.

        Looks for a MERlin output ``images`` subdirectory under the
        configured ``merlin_data_dir``.  If that does not exist, falls
        back to ``merlin_data_dir`` itself.
        """
        merlin_dir = Path(self.config.paths.merlin_data_dir)
        images_subdir = merlin_dir / "images"
        if images_subdir.is_dir():
            return images_subdir
        return merlin_dir

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
