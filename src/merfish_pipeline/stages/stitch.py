"""``stitch`` stage -- build tile mosaics from bead images.

For each imaging round (or z-slice, depending on configuration), this stage:

1. Reads the standardized position data produced by the ``index`` stage.
2. Locates bead-channel TIFF files in the raw data directory.
3. Groups images by imaging round or z-slice.
4. For each group, computes tile placement from stage positions (micron
   coordinates converted to pixel coordinates), applies optional per-
   microscope flips and transpose, and places every tile on a large canvas
   using max-projection for overlapping regions.
5. Writes each mosaic as a multi-page TIFF stack.

Grouping modes
--------------
``group_by = "ir"``
    One output file per imaging round, with z-slices stacked along the
    leading axis (``ZYX``).

``group_by = "z"``
    One output file per z-slice, with imaging rounds stacked along the
    leading axis (``TYX``).
"""

from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import tifffile

from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.io.tiff_io import read_tiff
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_filename(ir: int, fov: int, z: int) -> str:
    """Build the canonical merFISH tile filename (all indices are 1-based)."""
    return f"merFISH_{ir:02d}_{fov:03d}_{z:02d}.TIFF"


def _build_file_index(images_dir: Path) -> set[str]:
    """Build a set of TIFF filenames for O(1) membership checks."""
    return {
        f.name
        for f in images_dir.iterdir()
        if f.suffix.upper() == ".TIFF"
    }


def _detect_ir_count(file_index: set[str]) -> int:
    """Scan the file index to find the maximum imaging-round number."""
    pattern = re.compile(r"merFISH_(\d+)_\d+_\d+\.TIFF", re.IGNORECASE)
    max_ir = 0
    for name in file_index:
        m = pattern.match(name)
        if m:
            max_ir = max(max_ir, int(m.group(1)))
    return max_ir


def _detect_z_count(positions_df: "pd.DataFrame") -> int:  # noqa: F821
    """Count z-position columns in the normalised positions dataframe."""
    return sum(1 for col in positions_df.columns if col.startswith("z_position_"))


def _build_position_map(
    positions_df: "pd.DataFrame",  # noqa: F821
) -> dict[int, tuple[int, int]]:
    """Map tile_number -> (normalised_grid_x, normalised_grid_y).

    Grid coordinates are shifted so that the minimum values become 0,
    matching the reference algorithm.

    When ``grid_pos_x`` / ``grid_pos_y`` columns are present they are used
    directly (integer grid indices).  Otherwise the stage positions (micron
    floats) are quantized into a grid by estimating the tile pitch from the
    minimum nonzero gap between unique coordinate values.
    """
    pos_map: dict[int, tuple[int, int]] = {}

    if "grid_pos_x" in positions_df.columns:
        x_offset = int(positions_df["grid_pos_x"].min())
        y_offset = int(positions_df["grid_pos_y"].min())

        for _, row in positions_df.iterrows():
            tile = int(row["tile_number"])
            gx = int(row["grid_pos_x"]) - x_offset
            gy = int(row["grid_pos_y"]) - y_offset
            pos_map[tile] = (gx, gy)
    else:
        # Derive grid positions from stage coordinates by quantizing to
        # tile pitch (the minimum nonzero gap between unique positions).
        x_vals = positions_df["stage_pos_x"].values
        y_vals = positions_df["stage_pos_y"].values

        if len(x_vals) > 1:
            sorted_x = np.sort(np.unique(np.round(x_vals, 1)))
            sorted_y = np.sort(np.unique(np.round(y_vals, 1)))
            dx = np.diff(sorted_x)
            dy = np.diff(sorted_y)
            pitch_x = float(dx[dx > 0].min()) if len(dx[dx > 0]) else 1.0
            pitch_y = float(dy[dy > 0].min()) if len(dy[dy > 0]) else 1.0
        else:
            pitch_x = pitch_y = 1.0

        x_min = float(x_vals.min())
        y_min = float(y_vals.min())

        for _, row in positions_df.iterrows():
            tile = int(row["tile_number"])
            gx = round((float(row["stage_pos_x"]) - x_min) / pitch_x)
            gy = round((float(row["stage_pos_y"]) - y_min) / pitch_y)
            pos_map[tile] = (gx, gy)

    return pos_map


def _grid_extent(
    position_map: dict[int, tuple[int, int]],
) -> tuple[int, int]:
    """Return (grid_cols, grid_rows) from the position map."""
    if not position_map:
        return 0, 0
    max_x = max(gx for gx, _ in position_map.values())
    max_y = max(gy for _, gy in position_map.values())
    return max_x + 1, max_y + 1


def _load_and_transform(
    filepath: Path,
    flip_v: bool,
    flip_h: bool,
    transpose: bool,
) -> np.ndarray:
    """Read a single tile and apply microscope-specific transforms.

    The transform order matches the reference implementation:
    transpose first, then vertical flip, then horizontal flip.
    """
    tile = read_tiff(filepath)

    # Handle multi-page TIFFs -- take the first plane.
    if tile.ndim > 2:
        tile = tile[0]

    if transpose:
        tile = tile.T
    if flip_v:
        tile = tile[::-1, :]
    if flip_h:
        tile = tile[:, ::-1]

    # Ensure contiguous memory when any view operations were applied.
    if flip_v or flip_h or transpose:
        tile = np.ascontiguousarray(tile)

    return tile


def _detect_tile_dimensions(
    images_dir: Path,
    file_index: set[str],
    default_dims: list[int],
) -> tuple[int, int]:
    """Auto-detect tile (width, height) from a sample TIFF.

    Falls back to the microscope ``image_dimensions`` config if no file
    is available.
    """
    for name in file_index:
        sample_path = images_dir / name
        try:
            img = read_tiff(sample_path)
            if img.ndim > 2:
                img = img[0]
            h, w = img.shape
            return w, h
        except Exception:
            continue

    # Fallback
    return int(default_dims[0]), int(default_dims[1])


def _stitch_single_mosaic(
    images_dir: Path,
    ir: int,
    z: int,
    position_map: dict[int, tuple[int, int]],
    grid_cols: int,
    grid_rows: int,
    tile_w: int,
    tile_h: int,
    flip_v: bool,
    flip_h: bool,
    transpose: bool,
    file_index: set[str],
    num_workers: int = 8,
) -> np.ndarray:
    """Stitch all FOV tiles for one (IR, Z) pair into a single mosaic.

    Tiles are loaded in parallel using a thread pool.  Overlapping pixels
    are resolved via ``np.maximum`` (max-projection).
    """
    mosaic_w = grid_cols * tile_w
    mosaic_h = grid_rows * tile_h
    mosaic = np.zeros((mosaic_h, mosaic_w), dtype=np.uint16)

    # Identify which tiles actually exist on disk.
    tasks: list[tuple[int, int, int, Path]] = []
    for tile_num, (gx, gy) in position_map.items():
        fname = _build_filename(ir, tile_num, z)
        if fname in file_index:
            tasks.append((tile_num, gx, gy, images_dir / fname))

    if not tasks:
        return mosaic

    def _load(task: tuple[int, int, int, Path]) -> tuple[int, int, np.ndarray]:
        _, gx, gy, fpath = task
        tile = _load_and_transform(fpath, flip_v, flip_h, transpose)
        return gx, gy, tile

    effective_workers = min(num_workers, len(tasks))

    with ThreadPoolExecutor(max_workers=effective_workers) as pool:
        futures = [pool.submit(_load, t) for t in tasks]
        for future in as_completed(futures):
            gx, gy, tile = future.result()

            y0 = gy * tile_h
            x0 = gx * tile_w

            # Guard against unexpected tile sizes.
            actual_h, actual_w = tile.shape
            h = min(actual_h, tile_h, mosaic_h - y0)
            w = min(actual_w, tile_w, mosaic_w - x0)

            # Max-projection for overlapping regions.
            region = mosaic[y0 : y0 + h, x0 : x0 + w]
            np.maximum(region, tile[:h, :w], out=region)

    return mosaic


# ---------------------------------------------------------------------------
# Stage implementation
# ---------------------------------------------------------------------------


@register_stage("stitch")
class StitchStage(PipelineStage):
    """Build tile mosaics from bead-channel images."""

    description = "Stitch bead tiles into mosaics"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that all required inputs are present."""
        errors: list[str] = []

        # 1. Position data must exist (from index stage or explicit override).
        pos_path = self._resolve_position_file()
        if not pos_path.exists():
            errors.append(f"Position file not found: {pos_path}")

        # 2. Bead image directory must exist.
        images_dir = self._resolve_images_dir()
        if not images_dir.exists():
            errors.append(f"Bead images directory not found: {images_dir}")
        elif not images_dir.is_dir():
            errors.append(f"Bead images path is not a directory: {images_dir}")

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if at least one mosaic TIFF exists in the output dir."""
        out = self.get_output_dir()
        if not out.exists():
            return False
        tiffs = list(out.glob("*_mosaic.tiff"))
        return len(tiffs) > 0

    def run(self, dry_run: bool = False) -> StageResult:
        """Execute the stitch stage."""
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        output_dir.mkdir(parents=True, exist_ok=True)

        # ----------------------------------------------------------
        # 1. Read position data
        # ----------------------------------------------------------
        pos_path = self._resolve_position_file()
        self.logger.info("Reading position data from %s", pos_path)
        positions_df = read_sheet(pos_path)

        if positions_df.empty:
            return StageResult(
                status="failed",
                error=f"Position file is empty: {pos_path}",
            )

        position_map = _build_position_map(positions_df)
        grid_cols, grid_rows = _grid_extent(position_map)
        self.logger.info(
            "Position map: %d tiles, grid %d x %d",
            len(position_map), grid_cols, grid_rows,
        )

        # ----------------------------------------------------------
        # 2. Locate bead images and build file index
        # ----------------------------------------------------------
        images_dir = self._resolve_images_dir()
        self.logger.info("Scanning bead images in %s", images_dir)
        file_index = _build_file_index(images_dir)
        self.logger.info("Found %d TIFF files", len(file_index))

        if not file_index:
            return StageResult(
                status="failed",
                error=f"No TIFF files found in {images_dir}",
            )

        # ----------------------------------------------------------
        # 3. Detect experiment dimensions
        # ----------------------------------------------------------
        num_irs = _detect_ir_count(file_index)
        num_z = _detect_z_count(positions_df)
        if num_z == 0:
            # Fallback: scan filenames for max z index.
            z_pattern = re.compile(r"merFISH_\d+_\d+_(\d+)\.TIFF", re.IGNORECASE)
            num_z = max(
                (int(m.group(1)) for n in file_index if (m := z_pattern.match(n))),
                default=1,
            )

        tile_w, tile_h = _detect_tile_dimensions(
            images_dir, file_index, self.config.microscope.image_dimensions
        )
        self.logger.info(
            "Detected: %d IRs, %d z-slices, tile %d x %d px",
            num_irs, num_z, tile_w, tile_h,
        )

        # ----------------------------------------------------------
        # 4. Determine ranges and grouping
        # ----------------------------------------------------------
        stitch_cfg = self.config.stitch

        ir_range = stitch_cfg.ir_range
        if ir_range:
            ir_start, ir_end = ir_range[0], ir_range[1]
        else:
            ir_start, ir_end = 1, num_irs

        z_range = stitch_cfg.z_range
        if z_range:
            z_start, z_end = z_range[0], z_range[1]
        else:
            z_start, z_end = 1, num_z

        ir_list = list(range(ir_start, ir_end + 1))
        z_list = list(range(z_start, z_end + 1))

        group_by = stitch_cfg.group_by
        self.logger.info(
            "Grouping by %s | IRs %d-%d | Z %d-%d",
            group_by, ir_start, ir_end, z_start, z_end,
        )

        if dry_run:
            n_mosaics = len(ir_list) if group_by == "ir" else len(z_list)
            return StageResult(
                status="skipped",
                metadata={
                    "dry_run": True,
                    "group_by": group_by,
                    "ir_range": [ir_start, ir_end],
                    "z_range": [z_start, z_end],
                    "expected_mosaics": n_mosaics,
                    "tile_count": len(position_map),
                    "grid_size": [grid_cols, grid_rows],
                    "tile_size": [tile_w, tile_h],
                },
            )

        # Microscope transforms
        flip_h = self.config.microscope.flip_horizontal
        flip_v = self.config.microscope.flip_vertical
        do_transpose = self.config.microscope.transpose
        num_workers = self.config.execution.max_workers

        # ----------------------------------------------------------
        # 5. Build mosaics
        # ----------------------------------------------------------
        output_files: list[str] = []

        if group_by == "ir":
            output_files = self._generate_by_ir(
                images_dir=images_dir,
                ir_list=ir_list,
                z_list=z_list,
                position_map=position_map,
                grid_cols=grid_cols,
                grid_rows=grid_rows,
                tile_w=tile_w,
                tile_h=tile_h,
                flip_v=flip_v,
                flip_h=flip_h,
                transpose=do_transpose,
                file_index=file_index,
                output_dir=output_dir,
                num_workers=num_workers,
            )
        else:
            output_files = self._generate_by_z(
                images_dir=images_dir,
                ir_list=ir_list,
                z_list=z_list,
                position_map=position_map,
                grid_cols=grid_cols,
                grid_rows=grid_rows,
                tile_w=tile_w,
                tile_h=tile_h,
                flip_v=flip_v,
                flip_h=flip_h,
                transpose=do_transpose,
                file_index=file_index,
                output_dir=output_dir,
                num_workers=num_workers,
            )

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata={
                "group_by": group_by,
                "ir_range": [ir_start, ir_end],
                "z_range": [z_start, z_end],
                "mosaics_generated": len(output_files),
                "tile_count": len(position_map),
                "grid_size": [grid_cols, grid_rows],
                "tile_size": [tile_w, tile_h],
                "mosaic_size": [grid_cols * tile_w, grid_rows * tile_h],
            },
        )

        self.write_run_metadata(result, start_time)
        return result

    # ------------------------------------------------------------------
    # Internal: mosaic generation
    # ------------------------------------------------------------------

    def _generate_by_ir(
        self,
        images_dir: Path,
        ir_list: list[int],
        z_list: list[int],
        position_map: dict[int, tuple[int, int]],
        grid_cols: int,
        grid_rows: int,
        tile_w: int,
        tile_h: int,
        flip_v: bool,
        flip_h: bool,
        transpose: bool,
        file_index: set[str],
        output_dir: Path,
        num_workers: int,
    ) -> list[str]:
        """One multi-page TIFF per imaging round, z-slices stacked."""
        self.logger.info(
            "Generating mosaics grouped by IR (%d rounds, %d z-slices each)",
            len(ir_list), len(z_list),
        )
        output_files: list[str] = []

        for ir in ir_list:
            self.logger.info("Stitching IR %d ...", ir)
            stack: list[np.ndarray] = []

            for z in z_list:
                mosaic = _stitch_single_mosaic(
                    images_dir=images_dir,
                    ir=ir,
                    z=z,
                    position_map=position_map,
                    grid_cols=grid_cols,
                    grid_rows=grid_rows,
                    tile_w=tile_w,
                    tile_h=tile_h,
                    flip_v=flip_v,
                    flip_h=flip_h,
                    transpose=transpose,
                    file_index=file_index,
                    num_workers=num_workers,
                )
                stack.append(mosaic)

            stack_array = np.stack(stack, axis=0)
            out_path = output_dir / f"IR{ir:02d}_mosaic.tiff"
            tifffile.imwrite(
                str(out_path),
                stack_array,
                imagej=True,
                metadata={"axes": "ZYX"},
            )
            self.logger.info(
                "Saved %s (shape: %s)", out_path, stack_array.shape,
            )
            output_files.append(str(out_path))

        self.logger.info("Generated %d IR mosaic stacks", len(output_files))
        return output_files

    def _generate_by_z(
        self,
        images_dir: Path,
        ir_list: list[int],
        z_list: list[int],
        position_map: dict[int, tuple[int, int]],
        grid_cols: int,
        grid_rows: int,
        tile_w: int,
        tile_h: int,
        flip_v: bool,
        flip_h: bool,
        transpose: bool,
        file_index: set[str],
        output_dir: Path,
        num_workers: int,
    ) -> list[str]:
        """One multi-page TIFF per z-slice, imaging rounds stacked."""
        self.logger.info(
            "Generating mosaics grouped by Z (%d z-slices, %d rounds each)",
            len(z_list), len(ir_list),
        )
        output_files: list[str] = []

        for z in z_list:
            self.logger.info("Stitching Z %d ...", z)
            stack: list[np.ndarray] = []

            for ir in ir_list:
                mosaic = _stitch_single_mosaic(
                    images_dir=images_dir,
                    ir=ir,
                    z=z,
                    position_map=position_map,
                    grid_cols=grid_cols,
                    grid_rows=grid_rows,
                    tile_w=tile_w,
                    tile_h=tile_h,
                    flip_v=flip_v,
                    flip_h=flip_h,
                    transpose=transpose,
                    file_index=file_index,
                    num_workers=num_workers,
                )
                stack.append(mosaic)

            stack_array = np.stack(stack, axis=0)
            out_path = output_dir / f"Z{z:02d}_mosaic.tiff"
            tifffile.imwrite(
                str(out_path),
                stack_array,
                imagej=True,
                metadata={"axes": "TYX"},
            )
            self.logger.info(
                "Saved %s (shape: %s)", out_path, stack_array.shape,
            )
            output_files.append(str(out_path))

        self.logger.info("Generated %d Z-slice mosaic stacks", len(output_files))
        return output_files

    # ------------------------------------------------------------------
    # Internal: path resolution
    # ------------------------------------------------------------------

    def _resolve_position_file(self) -> Path:
        """Resolve the position file path.

        Priority:
        1. Explicit ``stitch.position_file`` config override.
        2. Standardized positions produced by the ``index`` stage.
        """
        if self.config.stitch.position_file is not None:
            return Path(self.config.stitch.position_file)
        return Path(self.config.paths.output_dir) / "index" / "positions.standardized.csv"

    def _resolve_images_dir(self) -> Path:
        """Resolve the bead images directory.

        Priority:
        1. Explicit ``stitch.images_dir`` config override.
        2. ``{raw_data_dir}/{bead_channel_folder}``
        """
        if self.config.stitch.images_dir is not None:
            return Path(self.config.stitch.images_dir)
        return Path(self.config.paths.raw_data_dir) / self.config.raw_data.bead_channel_folder
