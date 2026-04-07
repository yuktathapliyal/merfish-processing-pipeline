"""``filter_barcodes`` stage -- drop duplicated z-slices from MERlin barcode output.

When z-reregistration pads the z-stack to a uniform depth it can introduce
duplicated slices (tail-padding).  MERlin still decodes barcodes on those
planes, so this stage removes the corresponding rows from ``barcodes.csv``
to prevent inflated counts.

Algorithm
---------
1. Read ``barcodes.csv`` from MERlin output (contains FOV, z, barcode_id,
   mean_distance columns among others).
2. Read ``zmap_new_to_old.csv`` from the reregistration stage (maps new_z to
   old_z with an ``is_duplicate`` boolean flag).
3. Build per-FOV duplicate sets using the configured mode:

   - **"any"** (default): a ``new_z`` is duplicated if *any* IR has it
     flagged as duplicate (union across IRs).
   - **"all"**: a ``new_z`` is duplicated only if *all* IRs have it flagged
     as duplicate (intersection across IRs).

4. Convert 1-based ``new_z`` indices to 0-based for matching the barcodes
   table.
5. Filter barcode rows where ``(fov, z)`` matches a duplicate entry.
6. Write filtered barcodes CSV and removal reports.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

from merfish_pipeline.io.columns import FOV_CANDIDATES, Z_CANDIDATES
from merfish_pipeline.io.sheet_io import read_sheet, write_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)

# Required columns in the zmap CSV
_ZMAP_REQUIRED_COLS = {"FOV", "IR", "new_z", "old_z", "is_duplicate"}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _detect_columns(df: pd.DataFrame) -> tuple[str, str]:
    """Auto-detect the FOV and z column names in a barcodes DataFrame.

    Tries a list of common names (case-sensitive) and returns the first match
    for each.

    Returns
    -------
    (fov_col, z_col)

    Raises
    ------
    ValueError
        If either column cannot be detected.
    """
    fov_col: str | None = None
    z_col: str | None = None

    for candidate in FOV_CANDIDATES:
        if candidate in df.columns:
            fov_col = candidate
            break

    for candidate in Z_CANDIDATES:
        if candidate in df.columns:
            z_col = candidate
            break

    if fov_col is None or z_col is None:
        missing = []
        if fov_col is None:
            missing.append(f"FOV (tried {FOV_CANDIDATES})")
        if z_col is None:
            missing.append(f"z (tried {Z_CANDIDATES})")
        raise ValueError(
            f"Cannot auto-detect barcodes columns: {', '.join(missing)}. "
            f"Available columns: {list(df.columns)}"
        )

    return fov_col, z_col


def _build_duplicate_sets(
    zmap_df: pd.DataFrame, mode: str
) -> dict[int, set[int]]:
    """Compute per-FOV sets of duplicate z-indices (0-based).

    Parameters
    ----------
    zmap_df:
        DataFrame from ``zmap_new_to_old.csv`` with columns
        ``FOV, IR, new_z, old_z, is_duplicate``.
    mode:
        ``"any"`` -- union of duplicate ``new_z`` across IRs per FOV.
        ``"all"`` -- intersection of duplicate ``new_z`` across IRs per FOV.

    Returns
    -------
    dict mapping FOV (int) to a set of 0-based z-indices that are duplicated.
    """
    if not _ZMAP_REQUIRED_COLS.issubset(set(zmap_df.columns)):
        raise ValueError(
            "zmap_new_to_old.csv is missing required columns. "
            f"Expected at least: {sorted(_ZMAP_REQUIRED_COLS)}; "
            f"got: {list(zmap_df.columns)}"
        )

    zmap = zmap_df.copy()
    zmap["FOV"] = pd.to_numeric(zmap["FOV"], errors="raise").astype(int)
    zmap["IR"] = pd.to_numeric(zmap["IR"], errors="raise").astype(int)
    zmap["new_z"] = pd.to_numeric(zmap["new_z"], errors="raise").astype(int)
    zmap["is_duplicate"] = pd.to_numeric(
        zmap["is_duplicate"], errors="raise"
    ).astype(int)

    dup_sets: dict[int, set[int]] = {}

    for fov, fov_group in zmap.groupby("FOV"):
        ir_dups: dict[int, set[int]] = {
            ir: set(
                ir_group.loc[
                    ir_group["is_duplicate"] == 1, "new_z"
                ].tolist()
            )
            for ir, ir_group in fov_group.groupby("IR")
        }

        if not ir_dups:
            dup_sets[int(fov)] = set()
            continue

        if mode == "all":
            combined: set[int] | None = None
            for d in ir_dups.values():
                combined = d if combined is None else (combined & d)
            dup_1based = combined if combined is not None else set()
        else:
            # mode == "any" (default)
            dup_1based: set[int] = set()
            for d in ir_dups.values():
                dup_1based |= d

        # Convert 1-based new_z to 0-based for barcode matching
        dup_0based = {nz - 1 for nz in dup_1based if nz >= 1}
        dup_sets[int(fov)] = dup_0based

    return dup_sets


def _filter_barcodes(
    barcodes_df: pd.DataFrame,
    dup_sets: dict[int, set[int]],
    fov_col: str,
    z_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Filter barcode rows whose (FOV, z) falls in a duplicate set.

    Parameters
    ----------
    barcodes_df:
        The full barcodes DataFrame.
    dup_sets:
        Per-FOV sets of 0-based duplicate z-indices from
        :func:`_build_duplicate_sets`.
    fov_col:
        Name of the FOV column in the barcodes DataFrame.
    z_col:
        Name of the z column in the barcodes DataFrame.

    Returns
    -------
    (kept_df, removed_df)
        Two DataFrames: the rows that survive filtering and the rows that
        were removed.
    """
    df = barcodes_df.copy()
    df[fov_col] = pd.to_numeric(df[fov_col], errors="raise").astype(int)
    df[z_col] = pd.to_numeric(df[z_col], errors="raise").astype(int)

    is_dup = df.apply(
        lambda row: int(row[z_col]) in dup_sets.get(int(row[fov_col]), set()),
        axis=1,
    )

    kept_df = df.loc[~is_dup].copy()
    removed_df = df.loc[is_dup].copy()

    return kept_df, removed_df


def _write_removal_reports(
    removed_df: pd.DataFrame,
    output_dir: Path,
    fov_col: str,
    z_col: str,
    dup_sets: dict[int, set[int]],
    all_fovs: list[int],
    mode: str,
    zmap_df: pd.DataFrame,
) -> list[str]:
    """Write summary and detailed removal reports.

    Parameters
    ----------
    removed_df:
        Rows that were removed from the barcodes table.
    output_dir:
        Directory in which to write reports.
    fov_col:
        Name of the FOV column.
    z_col:
        Name of the z column.
    dup_sets:
        Per-FOV sets of 0-based duplicate z-indices.
    all_fovs:
        Sorted list of all FOV ids present in the original barcodes.
    mode:
        The duplicate aggregation mode that was used (``"any"`` or ``"all"``).
    zmap_df:
        The original zmap DataFrame for building the detailed report.

    Returns
    -------
    List of file paths written (as strings).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []

    # ---- Summary report: one row per FOV ----
    summary_rows = []
    for fov in all_fovs:
        z0_set = sorted(dup_sets.get(fov, set()))
        z1_set = [z + 1 for z in z0_set]
        summary_rows.append(
            {
                "FOV": fov,
                "removed_new_z0": ",".join(map(str, z0_set)),
                "removed_new_z1": ",".join(map(str, z1_set)),
                "count_removed": len(z0_set),
                "mode": mode,
            }
        )
    summary_path = output_dir / "removed_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    written.append(str(summary_path))

    # ---- Detailed report: per-(FOV, IR, new_z) for removed planes ----
    zmap_int = zmap_df.copy()
    for c in ["FOV", "IR", "new_z", "old_z", "is_duplicate"]:
        zmap_int[c] = pd.to_numeric(zmap_int[c], errors="raise").astype(int)
    zmap_int["new_z0"] = zmap_int["new_z"] - 1

    zmap_int["selected_for_removal"] = zmap_int.apply(
        lambda r: 1
        if int(r["new_z0"]) in dup_sets.get(int(r["FOV"]), set())
        else 0,
        axis=1,
    )

    detail_df = (
        zmap_int.loc[
            zmap_int["selected_for_removal"] == 1,
            ["FOV", "IR", "new_z0", "new_z", "old_z", "is_duplicate"],
        ]
        .rename(columns={"new_z": "new_z1"})
        .copy()
    )
    detail_df.insert(0, "mode", mode)
    detail_df = detail_df.sort_values(["FOV", "new_z0", "IR"])

    detail_path = output_dir / "removed_detailed.csv"
    detail_df.to_csv(detail_path, index=False)
    written.append(str(detail_path))

    return written


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("filter_barcodes")
class FilterBarcodesStage(PipelineStage):
    """Drop duplicated z-slices from MERlin barcode output after reregistration."""

    description = (
        "Filter MERlin barcodes to remove rows from duplicated z-slices "
        "introduced during z-reregistration padding"
    )

    # Expected output filenames
    _FILTERED_NAME = "barcodes_filtered.csv"
    _SUMMARY_NAME = "removed_summary.csv"
    _DETAILED_NAME = "removed_detailed.csv"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that the required input files exist."""
        errors: list[str] = []

        barcodes_path = self._resolve_barcodes_path()
        if barcodes_path is None:
            merlin_data_name = Path(self.config.paths.merlin_data_dir).name
            expected = (
                Path(self.config.paths.output_dir)
                / "merlin_analysis" / merlin_data_name / "ExportBarcodes" / "barcodes.csv"
            )
            errors.append(
                f"Cannot locate barcodes.csv. Expected at {expected}. "
                f"Set filter_barcodes.barcodes_file to the correct path "
                f"if MERlin output is elsewhere."
            )
        elif not barcodes_path.exists():
            errors.append(f"Barcodes file does not exist: {barcodes_path}")

        zmap_path = self._resolve_zmap_path()
        if zmap_path is None:
            errors.append(
                "Cannot locate zmap_new_to_old.csv. Provide "
                "filter_barcodes.zmap_file in the config or ensure "
                "reregistration output exists."
            )
        elif not zmap_path.exists():
            errors.append(f"Z-map file does not exist: {zmap_path}")

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if all expected output files already exist."""
        out = self.get_output_dir()
        return all(
            (out / name).exists()
            for name in (self._FILTERED_NAME, self._SUMMARY_NAME, self._DETAILED_NAME)
        )

    def run(self, dry_run: bool = False) -> StageResult:
        """Execute barcode filtering."""
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        output_dir.mkdir(parents=True, exist_ok=True)

        mode = self.config.filter_barcodes.mode

        barcodes_path = self._resolve_barcodes_path()
        zmap_path = self._resolve_zmap_path()

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would filter barcodes from %s using zmap %s (mode=%s)",
                barcodes_path,
                zmap_path,
                mode,
            )
            result = StageResult(
                status="skipped",
                metadata={"dry_run": True, "mode": mode},
            )
            self.write_run_metadata(result, start_time, parameters={"mode": mode})
            return result

        # ----------------------------------------------------------
        # Step 1: Read inputs
        # ----------------------------------------------------------
        self.logger.info("Reading barcodes from %s ...", barcodes_path)
        barcodes_df = read_sheet(barcodes_path)
        self.logger.info(
            "Loaded %d barcode rows with columns: %s",
            len(barcodes_df),
            list(barcodes_df.columns),
        )

        self.logger.info("Reading z-map from %s ...", zmap_path)
        zmap_df = read_sheet(zmap_path)
        self.logger.info(
            "Loaded %d z-map rows.", len(zmap_df)
        )

        # ----------------------------------------------------------
        # Step 2: Detect columns in barcodes
        # ----------------------------------------------------------
        try:
            fov_col, z_col = _detect_columns(barcodes_df)
        except ValueError as exc:
            return StageResult(
                status="failed",
                error=str(exc),
            )
        self.logger.info(
            "Detected barcode columns: fov=%r, z=%r", fov_col, z_col
        )

        # ----------------------------------------------------------
        # Step 3: Build per-FOV duplicate sets
        # ----------------------------------------------------------
        try:
            dup_sets = _build_duplicate_sets(zmap_df, mode=mode)
        except ValueError as exc:
            return StageResult(
                status="failed",
                error=str(exc),
            )

        total_dup_planes = sum(len(s) for s in dup_sets.values())
        self.logger.info(
            "Built duplicate sets for %d FOVs (%d total duplicate z-planes, mode=%s).",
            len(dup_sets),
            total_dup_planes,
            mode,
        )

        # ----------------------------------------------------------
        # Step 4: Filter barcodes
        # ----------------------------------------------------------
        kept_df, removed_df = _filter_barcodes(
            barcodes_df, dup_sets, fov_col, z_col
        )

        self.logger.info(
            "Filtering complete: %d rows kept, %d rows removed (of %d total).",
            len(kept_df),
            len(removed_df),
            len(barcodes_df),
        )

        # ----------------------------------------------------------
        # Step 5: Write outputs
        # ----------------------------------------------------------
        # Filtered barcodes
        filtered_path = output_dir / self._FILTERED_NAME
        write_sheet(kept_df, filtered_path)
        self.logger.info("Wrote filtered barcodes: %s", filtered_path)

        # All FOVs present in original barcodes (for complete reporting)
        barcodes_df[fov_col] = pd.to_numeric(
            barcodes_df[fov_col], errors="raise"
        ).astype(int)
        all_fovs = sorted(barcodes_df[fov_col].unique().tolist())

        # Removal reports
        report_files = _write_removal_reports(
            removed_df=removed_df,
            output_dir=output_dir,
            fov_col=fov_col,
            z_col=z_col,
            dup_sets=dup_sets,
            all_fovs=all_fovs,
            mode=mode,
            zmap_df=zmap_df,
        )
        for rp in report_files:
            self.logger.info("Wrote report: %s", rp)

        # ----------------------------------------------------------
        # Result
        # ----------------------------------------------------------
        output_files = [str(filtered_path)] + report_files

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata={
                "mode": mode,
                "rows_in": len(barcodes_df),
                "rows_kept": len(kept_df),
                "rows_removed": len(removed_df),
                "fovs_with_duplicates": sum(
                    1 for s in dup_sets.values() if s
                ),
                "total_duplicate_z_planes": total_dup_planes,
                "fov_col": fov_col,
                "z_col": z_col,
            },
        )

        self.write_run_metadata(
            result,
            start_time,
            parameters={"mode": mode},
        )

        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_barcodes_path(self) -> Path | None:
        """Resolve the path to the barcodes CSV.

        Priority:
        1. Explicit ``filter_barcodes.barcodes_file`` config override.
        2. Auto-detect at ``{output_dir}/merlin_analysis/{merlin_data_dir.name}/ExportBarcodes/barcodes.csv``.
        """
        explicit = self.config.filter_barcodes.barcodes_file
        if explicit is not None:
            return Path(explicit)

        # Auto-detect from MERlin analysis output
        output_dir = Path(self.config.paths.output_dir)
        merlin_data_name = Path(self.config.paths.merlin_data_dir).name
        candidate = (
            output_dir / "merlin_analysis" / merlin_data_name / "ExportBarcodes" / "barcodes.csv"
        )
        if candidate.exists():
            return candidate

        return None

    def _resolve_zmap_path(self) -> Path | None:
        """Resolve the path to the zmap_new_to_old CSV.

        Uses the explicit config value if provided; otherwise looks in the
        reregistration output directory.
        """
        explicit = self.config.filter_barcodes.zmap_file
        if explicit is not None:
            return Path(explicit)

        # Auto-detect from reregistration output
        output_dir = getattr(self.config.paths, "output_dir", None)
        if output_dir is not None:
            candidate = Path(output_dir) / "reregistration" / "zmap_new_to_old.csv"
            if candidate.exists():
                return candidate

        return None
