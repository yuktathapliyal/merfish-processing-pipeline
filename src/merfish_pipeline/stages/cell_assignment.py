"""``cell_assignment`` stage -- assign decoded barcodes to segmented cells.

This stage maps each barcode in the MERlin output to a cell by indexing into
per-FOV segmentation masks using the barcode's pixel coordinates.

Algorithm
---------
1. Load the barcodes CSV (from ``filter_barcodes``, explicit config, or MERlin
   output -- same fallback chain as the ``correlation`` stage).
2. Discover per-FOV mask TIFFs in the masks directory (from the
   ``segmentation`` stage output or an explicit override).
3. For each FOV present in the barcodes table:

   a. Load the matching mask TIFF.
   b. Round barcode ``x``, ``y``, ``z`` to integer and clip to valid mask
      bounds.
   c. **3-D masks** ``(Z, Y, X)``: ``label = mask[z, y, x]``.
      **2-D masks** ``(Y, X)``: ``label = mask[y, x]`` (z is ignored).
   d. Format the cell ID as ``Cell{fov}_{label}`` for non-zero labels
      (background label 0 maps to ``None``).

4. Optionally filter border cells: any cell that has at least one barcode
   within ``crop_margin`` pixels of a FOV edge is removed entirely.
5. Write output CSVs and a per-FOV summary.

Outputs
-------
- ``{output_dir}/cell_assignment/barcodes_assigned.csv``
- ``{output_dir}/cell_assignment/barcodes_assigned_filtered.csv`` (when
  ``crop_margin > 0``)
- ``{output_dir}/cell_assignment/assignment_summary.csv``
- ``{output_dir}/cell_assignment/run_metadata.json``
"""

from __future__ import annotations

import logging
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from merfish_pipeline.io.columns import (
    FOV_CANDIDATES,
    LOCAL_X_CANDIDATES,
    LOCAL_Y_CANDIDATES,
    Z_CANDIDATES,
    detect_column,
)
from merfish_pipeline.io.sheet_io import read_sheet, write_sheet
from merfish_pipeline.io.tiff_io import read_tiff
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# FOV-to-mask mapping
# ---------------------------------------------------------------------------

#: Regex to extract the last integer from a mask filename stem.
_FOV_NUMBER_RE = re.compile(r"(\d+)")


def _build_fov_mask_map(
    masks_dir: Path, pattern: str
) -> dict[int, Path]:
    """Build a mapping from FOV number to mask file path.

    FOV numbers are extracted from mask filenames by finding the *last*
    contiguous integer in the stem (e.g. ``aligned_images_003_masks.tif``
    yields FOV 3).
    """
    mask_map: dict[int, Path] = {}
    for p in sorted(masks_dir.glob(pattern)):
        matches = _FOV_NUMBER_RE.findall(p.stem)
        if not matches:
            logger.warning("Cannot extract FOV number from mask file: %s", p.name)
            continue
        # Use the last integer found (the one right before "_masks")
        fov = int(matches[-1])
        if fov in mask_map:
            logger.warning(
                "Duplicate mask for FOV %d: %s vs %s", fov, mask_map[fov].name, p.name
            )
        mask_map[fov] = p
    return mask_map


# ---------------------------------------------------------------------------
# Per-FOV assignment
# ---------------------------------------------------------------------------


def _assign_fov(
    df_fov: pd.DataFrame,
    mask: np.ndarray,
    fov: int,
    x_col: str,
    y_col: str,
    z_col: str,
    fov_id_format: str,
) -> pd.Series:
    """Assign cell IDs for all barcodes in one FOV.

    Rows with NaN coordinates (``x``/``y``, plus ``z`` for 3-D masks) are
    skipped -- ``.astype(int)`` raises on NaN values, and we don't want a
    handful of malformed rows to take down the whole stage.  Skipped rows
    receive ``None`` in the returned Series via index alignment in the
    caller.

    Parameters
    ----------
    df_fov:
        Barcodes belonging to a single FOV.
    mask:
        Segmentation mask -- 3-D ``(Z, Y, X)`` or 2-D ``(Y, X)``.
    fov:
        FOV number (used in the cell ID string).
    x_col, y_col, z_col:
        Column names for per-FOV pixel coordinates.
    fov_id_format:
        Format string with ``{fov}`` and ``{label}`` placeholders.

    Returns
    -------
    pd.Series
        Cell ID strings indexed by ``df_fov.index``.  Background pixels and
        rows with NaN coordinates map to ``None``.
    """
    # Drop rows with NaN coordinates BEFORE the int cast.  For 3-D masks
    # we additionally require a non-NaN z; for 2-D masks z is unused.
    coord_cols = [x_col, y_col]
    if mask.ndim == 3:
        coord_cols.append(z_col)
    valid = df_fov[coord_cols].notna().all(axis=1)
    n_dropped = int((~valid).sum())
    if n_dropped > 0:
        logger.warning(
            "FOV %d: dropping %d barcode row(s) with NaN coordinates.",
            fov,
            n_dropped,
        )
    df_valid = df_fov.loc[valid]

    if df_valid.empty:
        # Nothing to assign -- return an empty Series indexed by the
        # original df_fov so the caller's loc-assignment is a no-op.
        return pd.Series(
            pd.array([], dtype=object), index=df_valid.index
        )

    x_vals = df_valid[x_col].values.round().astype(int)
    y_vals = df_valid[y_col].values.round().astype(int)

    if mask.ndim == 3:
        z_vals = df_valid[z_col].values.round().astype(int)
        max_z, max_y, max_x = mask.shape
        z_vals = np.clip(z_vals, 0, max_z - 1)
        y_vals = np.clip(y_vals, 0, max_y - 1)
        x_vals = np.clip(x_vals, 0, max_x - 1)
        labels = mask[z_vals, y_vals, x_vals]
    elif mask.ndim == 2:
        max_y, max_x = mask.shape
        y_vals = np.clip(y_vals, 0, max_y - 1)
        x_vals = np.clip(x_vals, 0, max_x - 1)
        labels = mask[y_vals, x_vals]
    else:
        raise ValueError(f"Unexpected mask shape {mask.shape}; expected 2-D or 3-D")

    cell_ids = pd.array(
        [
            fov_id_format.format(fov=fov, label=int(lbl)) if lbl > 0 else None
            for lbl in labels
        ],
        dtype=object,
    )
    return pd.Series(cell_ids, index=df_valid.index)


# ---------------------------------------------------------------------------
# Border filtering
# ---------------------------------------------------------------------------


def _find_border_cells(
    df: pd.DataFrame,
    mask_shapes: dict[int, tuple[int, ...]],
    fov_col: str,
    x_col: str,
    y_col: str,
    cell_id_col: str,
    crop_margin: int,
) -> set[str]:
    """Identify cells with any barcode within *crop_margin* of a FOV edge.

    Parameters
    ----------
    df:
        Barcodes DataFrame with ``Cell_ID`` already assigned.
    mask_shapes:
        Mapping of FOV number to mask shape (used for width/height).
    fov_col, x_col, y_col, cell_id_col:
        Column names.
    crop_margin:
        Pixel distance from FOV edge.  Barcodes with ``x < margin``,
        ``y < margin``, ``x > (width - margin)``, or ``y > (height - margin)``
        are considered border barcodes.

    Returns
    -------
    Set of ``Cell_ID`` values to remove.
    """
    cells_to_remove: set[str] = set()

    for fov, shape in mask_shapes.items():
        # shape is (Z, Y, X) or (Y, X)
        h = shape[-2]
        w = shape[-1]

        fov_mask = df[fov_col] == fov
        fov_df = df.loc[fov_mask]
        if fov_df.empty:
            continue

        x_vals = fov_df[x_col].values
        y_vals = fov_df[y_col].values

        border = (
            (x_vals < crop_margin)
            | (y_vals < crop_margin)
            | (x_vals > w - crop_margin)
            | (y_vals > h - crop_margin)
        )

        border_ids = fov_df.loc[border, cell_id_col].dropna().unique()
        cells_to_remove.update(border_ids)

    return cells_to_remove


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("cell_assignment")
class CellAssignmentStage(PipelineStage):
    """Assign decoded barcodes to segmented cells using per-FOV masks."""

    description = "Assign barcodes to cells using segmentation masks"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []

        barcodes_path = self._resolve_barcodes_path()
        if barcodes_path is None:
            errors.append(
                "Cannot locate barcodes CSV. Set cell_assignment.barcodes_file "
                "or ensure filter_barcodes / MERlin output exists."
            )
        elif not barcodes_path.exists():
            errors.append(f"Barcodes file does not exist: {barcodes_path}")

        masks_dir = self._resolve_masks_dir()
        if masks_dir is None:
            errors.append(
                "Cannot locate masks directory. Set cell_assignment.masks_dir "
                "or run the segmentation stage first."
            )
        elif not masks_dir.is_dir():
            errors.append(f"Masks path is not a directory: {masks_dir}")
        else:
            pattern = self.config.cell_assignment.mask_pattern
            masks = list(masks_dir.glob(pattern))
            if not masks:
                errors.append(
                    f"No mask TIFFs found in {masks_dir} "
                    f"matching pattern '{pattern}'"
                )

        return errors

    def check_outputs_exist(self) -> bool:
        out = self.get_output_dir()
        return (out / "barcodes_assigned.csv").exists()

    def run(self, dry_run: bool = False) -> StageResult:
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        cfg = self.config.cell_assignment

        # ----------------------------------------------------------
        # 1. Resolve inputs
        # ----------------------------------------------------------
        barcodes_path = self._resolve_barcodes_path()
        masks_dir = self._resolve_masks_dir()

        if barcodes_path is None or masks_dir is None:
            return StageResult(
                status="failed",
                error="Could not resolve barcodes or masks path.",
            )

        mask_map = _build_fov_mask_map(masks_dir, cfg.mask_pattern)
        self.logger.info(
            "Found %d mask file(s) in %s", len(mask_map), masks_dir
        )

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would assign barcodes from %s using %d mask(s) "
                "from %s (crop_margin=%d)",
                barcodes_path,
                len(mask_map),
                masks_dir,
                cfg.crop_margin,
            )
            return StageResult(status="skipped", metadata={"dry_run": True})

        # ----------------------------------------------------------
        # 2. Load barcodes
        # ----------------------------------------------------------
        self.logger.info("Loading barcodes from %s ...", barcodes_path)
        barcodes_df = read_sheet(barcodes_path)
        self.logger.info(
            "Loaded %d barcode rows with columns: %s",
            len(barcodes_df),
            list(barcodes_df.columns),
        )

        # ----------------------------------------------------------
        # 3. Detect columns
        # ----------------------------------------------------------
        try:
            fov_col = detect_column(barcodes_df, FOV_CANDIDATES, "FOV")
            x_col = detect_column(barcodes_df, LOCAL_X_CANDIDATES, "x")
            y_col = detect_column(barcodes_df, LOCAL_Y_CANDIDATES, "y")
            z_col = detect_column(barcodes_df, Z_CANDIDATES, "z")
            barcodes_df[fov_col] = pd.to_numeric(
                barcodes_df[fov_col], errors="raise"
            ).astype(int)
        except (ValueError, TypeError) as exc:
            return StageResult(
                status="failed",
                error=f"Could not parse FOV column as integer: {exc}",
            )

        self.logger.info(
            "Detected columns: fov=%r, x=%r, y=%r, z=%r",
            fov_col,
            x_col,
            y_col,
            z_col,
        )

        # ----------------------------------------------------------
        # 4. Assign cell IDs per FOV
        # ----------------------------------------------------------
        cell_id_col = "Cell_ID"
        barcodes_df[cell_id_col] = None

        fovs = sorted(barcodes_df[fov_col].unique())
        mask_shapes: dict[int, tuple[int, ...]] = {}
        summary_rows: list[dict] = []
        n_assigned = 0
        n_no_mask = 0

        for fov in fovs:
            if fov not in mask_map:
                self.logger.warning("No mask found for FOV %d, skipping.", fov)
                n_no_mask += 1
                summary_rows.append(
                    {"fov": fov, "status": "no_mask", "n_barcodes": 0,
                     "n_assigned": 0, "n_cells": 0}
                )
                continue

            mask = read_tiff(mask_map[fov])
            mask_shapes[fov] = mask.shape

            fov_idx = barcodes_df[fov_col] == fov
            df_fov = barcodes_df.loc[fov_idx]

            cell_ids = _assign_fov(
                df_fov=df_fov,
                mask=mask,
                fov=fov,
                x_col=x_col,
                y_col=y_col,
                z_col=z_col,
                fov_id_format=cfg.fov_id_format,
            )

            barcodes_df.loc[fov_idx, cell_id_col] = cell_ids

            n_fov_assigned = cell_ids.notna().sum()
            n_fov_cells = cell_ids.dropna().nunique()
            n_assigned += n_fov_assigned

            summary_rows.append(
                {
                    "fov": fov,
                    "status": "ok",
                    "n_barcodes": len(df_fov),
                    "n_assigned": int(n_fov_assigned),
                    "n_cells": int(n_fov_cells),
                    "mask_shape": str(mask.shape),
                }
            )

            self.logger.debug(
                "FOV %d: %d/%d barcodes assigned to %d cells",
                fov,
                n_fov_assigned,
                len(df_fov),
                n_fov_cells,
            )

        self.logger.info(
            "Assignment complete: %d/%d barcodes assigned across %d FOVs "
            "(%d FOVs without masks).",
            n_assigned,
            len(barcodes_df),
            len(fovs),
            n_no_mask,
        )

        # ----------------------------------------------------------
        # 5. Write outputs
        # ----------------------------------------------------------
        output_dir.mkdir(parents=True, exist_ok=True)
        output_files: list[str] = []

        # Assigned barcodes
        assigned_path = output_dir / "barcodes_assigned.csv"
        write_sheet(barcodes_df, assigned_path)
        output_files.append(str(assigned_path))
        self.logger.info("Wrote assigned barcodes: %s", assigned_path)

        # Summary
        summary_df = pd.DataFrame(summary_rows)
        summary_path = output_dir / "assignment_summary.csv"
        write_sheet(summary_df, summary_path)
        output_files.append(str(summary_path))

        # ----------------------------------------------------------
        # 6. Border filtering (optional)
        # ----------------------------------------------------------
        n_border_cells = 0
        n_rows_removed = 0

        if cfg.crop_margin > 0 and mask_shapes:
            self.logger.info(
                "Filtering border cells (crop_margin=%d) ...", cfg.crop_margin
            )

            border_cells = _find_border_cells(
                df=barcodes_df,
                mask_shapes=mask_shapes,
                fov_col=fov_col,
                x_col=x_col,
                y_col=y_col,
                cell_id_col=cell_id_col,
                crop_margin=cfg.crop_margin,
            )
            n_border_cells = len(border_cells)

            if border_cells:
                filtered_df = barcodes_df[
                    ~barcodes_df[cell_id_col].isin(border_cells)
                ].copy()
                n_rows_removed = len(barcodes_df) - len(filtered_df)
            else:
                filtered_df = barcodes_df.copy()

            filtered_path = output_dir / "barcodes_assigned_filtered.csv"
            write_sheet(filtered_df, filtered_path)
            output_files.append(str(filtered_path))
            self.logger.info(
                "Border filtering: removed %d cells (%d rows), "
                "wrote %s",
                n_border_cells,
                n_rows_removed,
                filtered_path,
            )

        # ----------------------------------------------------------
        # 7. Result
        # ----------------------------------------------------------
        total_cells = barcodes_df[cell_id_col].dropna().nunique()

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata={
                "n_barcodes": len(barcodes_df),
                "n_assigned": n_assigned,
                "n_background": len(barcodes_df) - n_assigned,
                "n_cells_total": total_cells,
                "n_fovs": len(fovs),
                "n_fovs_no_mask": n_no_mask,
                "crop_margin": cfg.crop_margin,
                "n_border_cells_removed": n_border_cells,
                "n_rows_removed_by_border_filter": n_rows_removed,
            },
        )

        self.write_run_metadata(
            result,
            start_time,
            parameters={
                "barcodes_file": str(barcodes_path),
                "masks_dir": str(masks_dir),
                "mask_pattern": cfg.mask_pattern,
                "crop_margin": cfg.crop_margin,
                "fov_id_format": cfg.fov_id_format,
            },
        )

        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_barcodes_path(self) -> Path | None:
        """Resolve the barcodes CSV input.

        Priority:
        1. ``filter_barcodes`` stage output.
        2. Explicit ``cell_assignment.barcodes_file`` config.
        3. MERlin ``ExportBarcodes/barcodes.csv``.
        """
        output_dir = Path(self.config.paths.output_dir)

        # 1. filter_barcodes stage output
        filter_output = output_dir / "filter_barcodes" / "barcodes_filtered.csv"
        if filter_output.exists():
            return filter_output

        # 2. Explicit config
        explicit = self.config.cell_assignment.barcodes_file
        if explicit is not None:
            return Path(explicit)

        # 3. MERlin output
        merlin_data_name = Path(self.config.paths.merlin_data_dir).name
        merlin_output = (
            output_dir / "merlin_analysis" / merlin_data_name
            / "ExportBarcodes" / "barcodes.csv"
        )
        if merlin_output.exists():
            return merlin_output

        return None

    def _resolve_masks_dir(self) -> Path | None:
        """Resolve the masks directory.

        Priority:
        1. Explicit ``cell_assignment.masks_dir`` config.
        2. Auto-detect ``{output_dir}/segmentation/masks/``.
        """
        explicit = self.config.cell_assignment.masks_dir
        if explicit is not None:
            return Path(explicit)

        candidate = Path(self.config.paths.output_dir) / "segmentation" / "masks"
        if candidate.is_dir():
            self.logger.info("Auto-detected masks directory: %s", candidate)
            return candidate

        return None
