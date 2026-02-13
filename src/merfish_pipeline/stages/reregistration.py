"""``reregistration`` stage -- z-slice remapping for uniform depth across FOVs.

This stage wraps the z-remap algorithm from ``preprocessReRegistration_v2.py``.
It is an optional stage (controlled by ``config.reregistration.enabled``) that
remaps z-slices so every (FOV, IR) pair has a uniform number of z-planes,
starting from the best-focus slice identified by the ``focus_qc`` stage.

Algorithm
---------
1. Read ``best_focus_slices.csv`` (produced by ``focus_qc``) with columns
   ``FOV, IR01, IR02, ...`` where values are integer z-slice indices.
2. Detect IR columns from the CSV headers (pattern ``IR\\d+``).
3. For each (FOV, IR) pair the start z-slice is the best-focus z value.
4. Compute a uniform *target_z* depth:
   - For each FOV, find the maximum best-focus z across all IRs.
   - ``min_start`` = min of those per-FOV maxima across all FOVs.
   - ``target_z = total_z - min_start + 1`` (unless overridden in config).
5. Build pair plans: for each (FOV, IR), map ``new_z`` (1-based sequential) to
   ``old_z`` (starting from best-focus). If a pair runs out of real z-slices
   before reaching *target_z*, pad by duplicating the last available slice.
6. Execute file copies: for every channel directory, copy
   ``{channel}/merFISH_{IR}_{FOV}_{old_z}.TIFF`` to
   ``{remapped_data_dir}/{channel}/merFISH_{IR}_{FOV}_{new_z}.TIFF``.
"""

from __future__ import annotations

import csv
import re
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

# ---------------------------------------------------------------------------
# Filename / column patterns
# ---------------------------------------------------------------------------

#: Matches IR column names such as ``IR01``, ``IR02``, ``ir10``.
_IR_COL_RE = re.compile(r"^IR(\d+)$", re.IGNORECASE)

#: Matches the canonical merFISH filename pattern.
_TIFF_NAME_RE = re.compile(
    r"^merFISH_(\d{2})_(\d{3})_(\d{2})\.(?:tif|tiff)$", re.IGNORECASE
)


# ---------------------------------------------------------------------------
# Internal data structure for a single (FOV, IR) plan
# ---------------------------------------------------------------------------

@dataclass
class _PairPlan:
    """Z-remap plan for one (FOV, IR) pair."""

    fov: int
    ir: int
    start_old_z: int
    avail_z: int
    target_z: int
    dup_start_new_z: int | None
    zmap_new_to_old: list[int] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Module-level helper functions (extracted for testability)
# ---------------------------------------------------------------------------


def _detect_ir_columns(df: pd.DataFrame) -> list[tuple[int, str]]:
    """Find columns matching the ``IR\\d+`` pattern.

    Returns a sorted list of ``(ir_index, column_name)`` tuples.

    Raises
    ------
    ValueError
        If no IR columns are found.
    """
    result: list[tuple[int, str]] = []
    for col in df.columns:
        m = _IR_COL_RE.match(col)
        if m:
            result.append((int(m.group(1)), col))
    if not result:
        raise ValueError(
            "No IR columns found in CSV (expected names like IR01, IR02, ...)"
        )
    result.sort(key=lambda x: x[0])
    return result


def _infer_total_z(channel_dir: Path) -> int:
    """Scan filenames in *channel_dir* for the maximum z-slice index.

    Parameters
    ----------
    channel_dir:
        A single channel directory containing ``merFISH_*.TIFF`` files.

    Returns
    -------
    int
        The largest z-index found (1-based).

    Raises
    ------
    ValueError
        If no matching TIFF files are found.
    """
    max_z = 0
    for p in channel_dir.iterdir():
        if not p.is_file():
            continue
        m = _TIFF_NAME_RE.match(p.name)
        if m:
            z = int(m.group(3))
            if z > max_z:
                max_z = z
    if max_z == 0:
        raise ValueError(
            f"Could not infer total_z: no matching TIFF files in {channel_dir}"
        )
    return max_z


def _compute_target_z(
    df: pd.DataFrame,
    ir_cols: list[tuple[int, str]],
    total_z: int,
) -> tuple[int, int]:
    """Compute a uniform target z-depth from focus data.

    For each FOV the maximum best-focus z across all IRs gives the
    ``fov_start``.  The minimum ``fov_start`` across all FOVs gives
    ``min_start``.  The target depth is ``total_z - min_start + 1``.

    Parameters
    ----------
    df:
        Best-focus DataFrame with ``FOV`` and IR columns.
    ir_cols:
        Output of :func:`_detect_ir_columns`.
    total_z:
        Total number of z-slices in the raw data.

    Returns
    -------
    tuple[int, int]
        ``(target_z, min_start)``.

    Raises
    ------
    ValueError
        If the computed values are out of range.
    """
    col_names = [col for _, col in ir_cols]
    # Per-row (FOV) maximum of the IR values
    row_max = (
        pd.to_numeric(df[col_names].stack(), errors="raise")
        .unstack()
        .max(axis=1)
    )
    min_start = int(row_max.min())
    if min_start < 1:
        raise ValueError(
            f"Computed min_start < 1 from focus CSV (got {min_start})"
        )
    target_z = total_z - min_start + 1
    if target_z < 1:
        raise ValueError(
            f"Computed target_z < 1 (total_z={total_z}, min_start={min_start})"
        )
    return target_z, min_start


def _build_pair_plans(
    df: pd.DataFrame,
    ir_cols: list[tuple[int, str]],
    total_z: int,
    target_z: int,
) -> dict[tuple[int, int], _PairPlan]:
    """Build a z-remap plan for every (FOV, IR) pair.

    Parameters
    ----------
    df:
        Best-focus DataFrame with ``FOV`` column and IR columns.
    ir_cols:
        Sorted list of ``(ir_index, column_name)`` from :func:`_detect_ir_columns`.
    total_z:
        Total number of original z-slices per stack.
    target_z:
        Desired uniform depth for all output stacks.

    Returns
    -------
    dict[tuple[int, int], _PairPlan]
        Keyed by ``(fov, ir_index)``.

    Raises
    ------
    ValueError
        On non-numeric or out-of-range best-focus values.
    """
    if "FOV" not in df.columns:
        raise ValueError("best_focus_slices.csv must contain a 'FOV' column")

    plans: dict[tuple[int, int], _PairPlan] = {}

    for _, row in df.iterrows():
        fov_raw = row["FOV"]
        try:
            fov = int(fov_raw)
        except Exception:
            raise ValueError(f"FOV value '{fov_raw}' is not an integer")

        for ir_idx, col_name in ir_cols:
            try:
                start_old_z = int(pd.to_numeric(row[col_name], errors="raise"))
            except Exception as exc:
                raise ValueError(
                    f"Non-numeric best-focus value for FOV {fov}, {col_name}: {exc}"
                )

            if start_old_z < 1 or start_old_z > total_z:
                raise ValueError(
                    f"start_old_z out of range for FOV {fov}, {col_name}: "
                    f"{start_old_z} (expected 1..{total_z})"
                )

            avail_z = max(0, total_z - start_old_z + 1)

            # Build new_z -> old_z mapping (length == target_z)
            last_real_old = min(start_old_z + max(0, avail_z - 1), total_z)
            zmap: list[int] = []
            for new_z in range(1, target_z + 1):
                old_z = start_old_z + (new_z - 1)
                if old_z > total_z:
                    old_z = total_z
                if new_z > avail_z:
                    old_z = last_real_old
                zmap.append(int(old_z))

            dup_start = (avail_z + 1) if avail_z < target_z else None

            plans[(fov, ir_idx)] = _PairPlan(
                fov=fov,
                ir=ir_idx,
                start_old_z=start_old_z,
                avail_z=avail_z,
                target_z=target_z,
                dup_start_new_z=int(dup_start) if dup_start is not None else None,
                zmap_new_to_old=zmap,
            )

    return plans


def _resolve_src_path(
    ch_dir: Path, ir: int, fov: int, z: int
) -> Path | None:
    """Find the source TIFF, tolerating ``.TIFF`` / ``.tif`` extensions."""
    base = f"merFISH_{ir:02d}_{fov:03d}_{z:02d}"
    for ext in (".TIFF", ".tif"):
        candidate = ch_dir / f"{base}{ext}"
        if candidate.exists() and candidate.is_file():
            return candidate
    return None


def _execute_copies(
    plans: dict[tuple[int, int], _PairPlan],
    channels: list[str],
    in_dir: Path,
    out_dir: Path,
    logger: Any,
    strict_missing: bool = False,
    dry_run: bool = False,
) -> tuple[int, int, int]:
    """Perform (or simulate) all file copies described by *plans*.

    Parameters
    ----------
    plans:
        Per-(FOV, IR) remap plans from :func:`_build_pair_plans`.
    channels:
        Channel folder names (e.g. ``["488nm, Raw", "561nm, Raw"]``).
    in_dir:
        Root input directory (``raw_data_dir``).
    out_dir:
        Root output directory (``remapped_data_dir``).
    logger:
        Logger instance for progress messages.
    strict_missing:
        If ``True``, raise :class:`FileNotFoundError` on a missing source file
        instead of skipping it.
    dry_run:
        If ``True``, do not perform any actual file I/O.

    Returns
    -------
    tuple[int, int, int]
        ``(n_copied, n_skipped, n_total)`` counts.
    """
    # Ensure output channel directories exist
    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        for ch in channels:
            (out_dir / ch).mkdir(parents=True, exist_ok=True)

    n_copied = 0
    n_skipped = 0
    n_total = 0

    sorted_keys = sorted(plans.keys())
    total_pairs = len(sorted_keys)

    for idx, (fov, ir) in enumerate(sorted_keys, 1):
        plan = plans[(fov, ir)]

        for new_z, old_z in enumerate(plan.zmap_new_to_old, start=1):
            for ch in channels:
                n_total += 1
                ch_dir = in_dir / ch
                src = _resolve_src_path(ch_dir, ir, fov, old_z)

                dst = (
                    out_dir
                    / ch
                    / f"merFISH_{ir:02d}_{fov:03d}_{new_z:02d}.TIFF"
                )

                if src is None:
                    msg = (
                        f"Missing source: "
                        f"{ch_dir / f'merFISH_{ir:02d}_{fov:03d}_{old_z:02d}.TIFF'}"
                    )
                    if strict_missing:
                        raise FileNotFoundError(msg)
                    logger.warning(msg)
                    n_skipped += 1
                    continue

                if not dry_run:
                    shutil.copy2(str(src), str(dst))
                n_copied += 1

        if idx % 50 == 0 or idx == total_pairs:
            logger.info(
                "  Copy progress: %d / %d (FOV, IR) pairs ...", idx, total_pairs
            )

    return n_copied, n_skipped, n_total


# ---------------------------------------------------------------------------
# CSV / metadata writers
# ---------------------------------------------------------------------------


def _write_zmap_csv(
    plans: dict[tuple[int, int], _PairPlan], path: Path
) -> None:
    """Write the long-format new_z -> old_z mapping CSV."""
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["FOV", "IR", "new_z", "old_z", "is_duplicate"])
        for fov, ir in sorted(plans.keys()):
            plan = plans[(fov, ir)]
            for new_z, old_z in enumerate(plan.zmap_new_to_old, start=1):
                is_dup = 1 if new_z > plan.avail_z else 0
                writer.writerow(
                    [
                        f"{plan.fov}",
                        f"{plan.ir:02d}",
                        f"{new_z:02d}",
                        f"{old_z:02d}",
                        is_dup,
                    ]
                )


def _write_per_fov_summary(
    df: pd.DataFrame,
    ir_cols: list[tuple[int, str]],
    total_z: int,
    target_z: int,
    path: Path,
) -> None:
    """Write per-FOV summary CSV with fov_start and limiting flag."""
    col_names = [col for _, col in ir_cols]
    row_max = (
        pd.to_numeric(df[col_names].stack(), errors="raise")
        .unstack()
        .max(axis=1)
    )
    per_fov = pd.DataFrame(
        {
            "FOV": pd.to_numeric(df["FOV"], errors="raise").astype(int),
            "fov_start": row_max.astype(int),
        }
    )
    per_fov["avail_z_at_fov_start"] = (total_z - per_fov["fov_start"] + 1).astype(int)
    min_start = int(per_fov["fov_start"].min())
    per_fov["target_z"] = int(target_z)
    per_fov["is_limiting"] = (per_fov["fov_start"] == min_start).astype(int)
    per_fov[
        ["FOV", "fov_start", "avail_z_at_fov_start", "target_z", "is_limiting"]
    ].sort_values("FOV").to_csv(path, index=False)


def _write_planned_copies(
    plans: dict[tuple[int, int], _PairPlan],
    channels: list[str],
    in_dir: Path,
    out_dir: Path,
    path: Path,
) -> None:
    """Write a CSV listing every planned copy operation."""
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "channel",
                "IR",
                "FOV",
                "new_z",
                "old_z",
                "src_path",
                "dst_path",
                "action",
            ]
        )
        for fov, ir in sorted(plans.keys()):
            plan = plans[(fov, ir)]
            for new_z, old_z in enumerate(plan.zmap_new_to_old, start=1):
                action = "COPY" if new_z <= plan.avail_z else "DUPLICATE"
                for ch in channels:
                    ch_dir = in_dir / ch
                    resolved = _resolve_src_path(ch_dir, ir, fov, old_z)
                    src_path = (
                        resolved
                        if resolved is not None
                        else ch_dir / f"merFISH_{ir:02d}_{fov:03d}_{old_z:02d}.TIFF"
                    )
                    dst_path = (
                        out_dir
                        / ch
                        / f"merFISH_{ir:02d}_{fov:03d}_{new_z:02d}.TIFF"
                    )
                    writer.writerow(
                        [
                            ch,
                            f"{ir:02d}",
                            f"{fov:03d}",
                            f"{new_z:02d}",
                            f"{old_z:02d}",
                            str(src_path),
                            str(dst_path),
                            action,
                        ]
                    )


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("reregistration")
class ReregistrationStage(PipelineStage):
    """Z-slice remapping for uniform depth across FOVs."""

    description = "Remap z-slices to create uniform depth across FOVs"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that required inputs exist."""
        errors: list[str] = []

        # The stage depends on focus_qc output
        focus_csv = self._focus_csv_path()
        if not focus_csv.exists():
            errors.append(
                f"best_focus_slices.csv not found at {focus_csv}. "
                f"Run the focus_qc stage first."
            )

        raw_dir = Path(self.config.paths.raw_data_dir)
        if not raw_dir.exists():
            errors.append(f"Raw data directory does not exist: {raw_dir}")
        elif not raw_dir.is_dir():
            errors.append(f"Raw data path is not a directory: {raw_dir}")
        else:
            channels = self._detect_channels()
            if not channels:
                errors.append(
                    f"No channel directories (containing 'nm') found in {raw_dir}"
                )

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if the run metadata and remapped data directory exist."""
        stage_dir = self.get_output_dir()
        metadata_file = stage_dir / "run_metadata.json"
        zmap_file = stage_dir / "zmap_new_to_old.csv"
        return metadata_file.exists() and zmap_file.exists()

    def run(self, dry_run: bool = False) -> StageResult:
        """Execute the z-remap algorithm."""
        start_time = datetime.now()
        stage_dir = self.get_output_dir()
        stage_dir.mkdir(parents=True, exist_ok=True)

        rereg_cfg = self.config.reregistration

        # Merge dry_run from config and method argument
        effective_dry_run = dry_run or rereg_cfg.dry_run

        # ----------------------------------------------------------
        # Step 1: Read best-focus CSV
        # ----------------------------------------------------------
        focus_csv = self._focus_csv_path()
        self.logger.info("Reading focus data from %s", focus_csv)
        df = pd.read_csv(focus_csv)

        # ----------------------------------------------------------
        # Step 2: Detect IR columns
        # ----------------------------------------------------------
        ir_cols = _detect_ir_columns(df)
        ir_indices = [idx for idx, _ in ir_cols]
        self.logger.info("Detected IR indices: %s", ir_indices)

        # ----------------------------------------------------------
        # Step 3: Auto-detect channels and total_z
        # ----------------------------------------------------------
        channels = self._detect_channels()
        if not channels:
            return StageResult(
                status="failed",
                error=(
                    f"No channel directories (containing 'nm') found "
                    f"in {self.config.paths.raw_data_dir}"
                ),
            )
        self.logger.info("Channels: %s", channels)

        raw_dir = Path(self.config.paths.raw_data_dir)

        if rereg_cfg.total_z is not None:
            total_z = rereg_cfg.total_z
            self.logger.info("total_z (config override): %d", total_z)
        else:
            # Infer from the first channel directory
            first_ch_dir = raw_dir / channels[0]
            total_z = _infer_total_z(first_ch_dir)
            self.logger.info("total_z (auto-detected): %d", total_z)

        # ----------------------------------------------------------
        # Step 4: Compute target_z
        # ----------------------------------------------------------
        if rereg_cfg.target_z is not None:
            target_z = rereg_cfg.target_z
            _, min_start = _compute_target_z(df, ir_cols, total_z)
            self.logger.info(
                "target_z (config override): %d  (min_start=%d)",
                target_z,
                min_start,
            )
        else:
            target_z, min_start = _compute_target_z(df, ir_cols, total_z)
            self.logger.info(
                "target_z (auto-computed): %d  (min_start=%d)",
                target_z,
                min_start,
            )

        if target_z < 1:
            return StageResult(
                status="failed",
                error=f"target_z must be >= 1 (got {target_z})",
            )

        # ----------------------------------------------------------
        # Step 5: Build pair plans
        # ----------------------------------------------------------
        self.logger.info("Building z-remap plans ...")
        plans = _build_pair_plans(df, ir_cols, total_z, target_z)
        self.logger.info("Built plans for %d (FOV, IR) pairs.", len(plans))

        # ----------------------------------------------------------
        # Step 6: Write diagnostic CSVs (always, even in dry-run)
        # ----------------------------------------------------------
        zmap_path = stage_dir / "zmap_new_to_old.csv"
        _write_zmap_csv(plans, zmap_path)
        self.logger.info("Wrote zmap: %s", zmap_path)

        per_fov_path = stage_dir / "per_fov_summary.csv"
        _write_per_fov_summary(df, ir_cols, total_z, target_z, per_fov_path)
        self.logger.info("Wrote per-FOV summary: %s", per_fov_path)

        out_dir = Path(self.config.paths.remapped_data_dir)
        planned_path = stage_dir / "planned_copies.csv"
        _write_planned_copies(plans, channels, raw_dir, out_dir, planned_path)
        self.logger.info("Wrote planned copies: %s", planned_path)

        output_files = [
            str(zmap_path),
            str(per_fov_path),
            str(planned_path),
        ]

        # ----------------------------------------------------------
        # Step 7: Execute file copies (unless dry-run)
        # ----------------------------------------------------------
        if effective_dry_run:
            self.logger.info(
                "[DRY RUN] Skipping file copies. %d planned operations.",
                sum(len(p.zmap_new_to_old) * len(channels) for p in plans.values()),
            )
            result = StageResult(
                status="completed",
                output_files=output_files,
                metadata={
                    "dry_run": True,
                    "total_z": total_z,
                    "target_z": target_z,
                    "min_start": min_start,
                    "n_pairs": len(plans),
                    "n_channels": len(channels),
                    "channels": channels,
                    "ir_indices": ir_indices,
                },
            )
            self.write_run_metadata(
                result,
                start_time,
                parameters={
                    "total_z": total_z,
                    "target_z": target_z,
                    "min_start": min_start,
                    "channels": channels,
                    "ir_indices": ir_indices,
                    "dry_run": True,
                },
            )
            return result

        self.logger.info("Executing file copies to %s ...", out_dir)
        n_copied, n_skipped, n_total = _execute_copies(
            plans=plans,
            channels=channels,
            in_dir=raw_dir,
            out_dir=out_dir,
            logger=self.logger,
            strict_missing=rereg_cfg.strict_missing,
            dry_run=False,
        )
        self.logger.info(
            "Copy complete: %d copied, %d skipped, %d total.",
            n_copied,
            n_skipped,
            n_total,
        )

        # ----------------------------------------------------------
        # Result
        # ----------------------------------------------------------
        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata={
                "dry_run": False,
                "total_z": total_z,
                "target_z": target_z,
                "min_start": min_start,
                "n_pairs": len(plans),
                "n_channels": len(channels),
                "channels": channels,
                "ir_indices": ir_indices,
                "n_copied": n_copied,
                "n_skipped": n_skipped,
                "n_total": n_total,
            },
        )

        self.write_run_metadata(
            result,
            start_time,
            parameters={
                "total_z": total_z,
                "target_z": target_z,
                "min_start": min_start,
                "channels": channels,
                "ir_indices": ir_indices,
                "strict_missing": rereg_cfg.strict_missing,
            },
        )

        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _focus_csv_path(self) -> Path:
        """Return the path to the best-focus CSV produced by the focus_qc stage."""
        return Path(self.config.paths.output_dir) / "focus_qc" / "best_focus_slices.csv"

    def _detect_channels(self) -> list[str]:
        """Auto-detect channel directories by scanning raw_data_dir for folders containing 'nm'."""
        raw_dir = Path(self.config.paths.raw_data_dir)
        if not raw_dir.is_dir():
            return []
        return sorted(
            p.name
            for p in raw_dir.iterdir()
            if p.is_dir() and "nm" in p.name.lower()
        )
