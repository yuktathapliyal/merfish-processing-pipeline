"""``inspect_positions`` stage -- analyse position drift and microscope logs.

Reads the standardized position data produced by the ``index`` stage and
computes per-round, per-FOV drift statistics.  Optionally parses the
microscope log file for laser power, exposure, and z-stack parameters.

Outputs
-------
drift_report.csv
    One row per (round, fov) pair with columns:
    ``round, fov, delta_x, delta_y, displacement``
drift_summary.txt
    Aggregate statistics (mean / std / max drift per round, plus any
    microscope-log parameters that were parsed).
drift_plot.png  *(optional)*
    Scatter plot of FOV positions coloured by round.  Only generated when
    *matplotlib* is importable.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal data structures
# ---------------------------------------------------------------------------

@dataclass
class _ImagingParams:
    """Parameters parsed from a microscope log file (e.g. ``merfish_log.txt``)."""

    z_stacks: int = 0
    upper_planes: int = 0
    lower_planes: int = 0
    initial_focus: float = 0.0

    laser_488_power: float = 0.0
    laser_561_power: float = 0.0
    laser_640_power: float = 0.0
    laser_750_power: float = 0.0

    exposure_488: float = 0.0
    exposure_561: float = 0.0
    exposure_640: float = 0.0
    exposure_750: float = 0.0


# ---------------------------------------------------------------------------
# Pure-function helpers (no side-effects, easily testable)
# ---------------------------------------------------------------------------

def _parse_microscope_log(log_path: Path) -> _ImagingParams:
    """Parse the header of a microscope log file for imaging parameters.

    The parser reads the first 30 lines (parameters are typically in the first
    ~15) and extracts laser powers, exposure times, z-stack counts, and focus
    settings.

    Returns a populated :class:`_ImagingParams` instance.  If the file cannot
    be read or a line cannot be parsed, defaults (zeros) are kept for those
    fields.
    """
    params = _ImagingParams()

    if not log_path.exists():
        return params

    try:
        text = log_path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()[:30]
    except Exception:
        return params

    _SIMPLE_INT = {
        "# Z stacks:": "z_stacks",
        "Upper Planes:": "upper_planes",
        "Lower Planes:": "lower_planes",
    }

    _SIMPLE_FLOAT = {
        "Initial focus position:": "initial_focus",
        "488nm laser power:": "laser_488_power",
        "561nm laser power:": "laser_561_power",
        "640nm laser power:": "laser_640_power",
        "750nm laser power:": "laser_750_power",
    }

    _EXPOSURE_FLOAT = {
        "488nm exposure:": "exposure_488",
        "561nm exposure:": "exposure_561",
        "640nm exposure:": "exposure_640",
        "750nm exposure:": "exposure_750",
    }

    for line in lines:
        stripped = line.strip()

        for prefix, attr in _SIMPLE_INT.items():
            if stripped.startswith(prefix):
                try:
                    setattr(params, attr, int(stripped.split(":", 1)[1].strip()))
                except (ValueError, IndexError):
                    pass
                break

        for prefix, attr in _SIMPLE_FLOAT.items():
            if stripped.startswith(prefix):
                try:
                    setattr(params, attr, float(stripped.split(":", 1)[1].strip()))
                except (ValueError, IndexError):
                    pass
                break

        for prefix, attr in _EXPOSURE_FLOAT.items():
            if stripped.startswith(prefix):
                try:
                    raw = stripped.split(":", 1)[1]
                    raw = raw.replace("ms", "").strip()
                    setattr(params, attr, float(raw))
                except (ValueError, IndexError):
                    pass
                break

    return params


def _compute_drift_report(
    positions_df: pd.DataFrame,
    rounds_to_check: Optional[list[int]] = None,
) -> pd.DataFrame:
    """Compute per-FOV displacement between each round and the reference round.

    The first round (numerically smallest) serves as the reference.  For every
    subsequent round the Euclidean displacement of each FOV (identified by
    ``tile_number``) is computed relative to that reference.

    Parameters
    ----------
    positions_df:
        Normalised position data with at least the columns
        ``round, tile_number, stage_pos_x, stage_pos_y``.
    rounds_to_check:
        Optional subset of round numbers to include.  When *None* all rounds
        in the data are used.

    Returns
    -------
    pd.DataFrame
        Columns: ``round, fov, delta_x, delta_y, displacement``.
        The reference round itself is included with zero-displacement rows so
        that every FOV in every round appears.
    """
    required = {"round", "tile_number", "stage_pos_x", "stage_pos_y"}
    if not required.issubset(positions_df.columns):
        missing = required - set(positions_df.columns)
        raise ValueError(
            f"Position data is missing required columns: {missing}"
        )

    df = positions_df.copy()

    # Optional round filtering
    if rounds_to_check is not None:
        df = df[df["round"].isin(rounds_to_check)]

    if df.empty:
        return pd.DataFrame(columns=["round", "fov", "delta_x", "delta_y", "displacement"])

    all_rounds = sorted(df["round"].unique())
    ref_round = all_rounds[0]
    ref = df[df["round"] == ref_round][["tile_number", "stage_pos_x", "stage_pos_y"]].copy()
    ref = ref.rename(columns={"stage_pos_x": "ref_x", "stage_pos_y": "ref_y"})

    records: list[dict[str, Any]] = []

    for rnd in all_rounds:
        rnd_df = df[df["round"] == rnd][["tile_number", "stage_pos_x", "stage_pos_y"]]
        merged = rnd_df.merge(ref, on="tile_number", how="inner")

        delta_x = merged["stage_pos_x"] - merged["ref_x"]
        delta_y = merged["stage_pos_y"] - merged["ref_y"]
        displacement = np.sqrt(delta_x ** 2 + delta_y ** 2)

        for tile, dx, dy, disp in zip(
            merged["tile_number"], delta_x, delta_y, displacement
        ):
            records.append(
                {
                    "round": rnd,
                    "fov": tile,
                    "delta_x": float(dx),
                    "delta_y": float(dy),
                    "displacement": float(disp),
                }
            )

    return pd.DataFrame(records, columns=["round", "fov", "delta_x", "delta_y", "displacement"])


def _format_drift_summary(
    drift_df: pd.DataFrame,
    imaging_params: Optional[_ImagingParams] = None,
) -> str:
    """Build a human-readable drift summary from the drift report.

    Parameters
    ----------
    drift_df:
        The output of :func:`_compute_drift_report`.
    imaging_params:
        Optional parsed microscope parameters to append at the bottom.

    Returns
    -------
    str
        Multi-line text suitable for writing to ``drift_summary.txt``.
    """
    lines: list[str] = []
    lines.append("MERFISH Position Drift Summary")
    lines.append("=" * 50)
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    if drift_df.empty:
        lines.append("No drift data available (empty position table).")
        return "\n".join(lines)

    # -- Overall statistics ---------------------------------------------------
    lines.append("Overall Drift Statistics")
    lines.append("-" * 50)
    lines.append(f"  Total FOV-round pairs : {len(drift_df)}")
    lines.append(f"  Rounds analysed       : {sorted(drift_df['round'].unique())}")
    lines.append(f"  Unique FOVs           : {drift_df['fov'].nunique()}")
    lines.append("")
    lines.append(f"  Mean displacement     : {drift_df['displacement'].mean():.6f}")
    lines.append(f"  Std displacement      : {drift_df['displacement'].std():.6f}")
    lines.append(f"  Max displacement      : {drift_df['displacement'].max():.6f}")
    lines.append(f"  Mean |delta_x|        : {drift_df['delta_x'].abs().mean():.6f}")
    lines.append(f"  Mean |delta_y|        : {drift_df['delta_y'].abs().mean():.6f}")
    lines.append("")

    # -- Per-round breakdown --------------------------------------------------
    lines.append("Per-Round Drift Statistics")
    lines.append("-" * 50)
    lines.append(
        f"{'Round':>6}  {'N_FOVs':>6}  {'Mean Disp':>10}  {'Std Disp':>10}  "
        f"{'Max Disp':>10}  {'Max |dX|':>10}  {'Max |dY|':>10}"
    )

    for rnd, grp in drift_df.groupby("round"):
        lines.append(
            f"{rnd:>6}  {len(grp):>6}  {grp['displacement'].mean():>10.6f}  "
            f"{grp['displacement'].std():>10.6f}  {grp['displacement'].max():>10.6f}  "
            f"{grp['delta_x'].abs().max():>10.6f}  {grp['delta_y'].abs().max():>10.6f}"
        )

    lines.append("")

    # -- Microscope log parameters (optional) ---------------------------------
    if imaging_params is not None and (
        imaging_params.z_stacks > 0
        or imaging_params.laser_488_power > 0
        or imaging_params.exposure_488 > 0
    ):
        lines.append("Microscope Log Parameters")
        lines.append("-" * 50)
        lines.append(f"  Z-stacks          : {imaging_params.z_stacks}")
        lines.append(f"  Upper planes      : {imaging_params.upper_planes}")
        lines.append(f"  Lower planes      : {imaging_params.lower_planes}")
        lines.append(f"  Initial focus     : {imaging_params.initial_focus}")
        lines.append("")
        lines.append("  Laser Powers (W):")
        lines.append(f"    488nm : {imaging_params.laser_488_power}")
        lines.append(f"    561nm : {imaging_params.laser_561_power}")
        lines.append(f"    640nm : {imaging_params.laser_640_power}")
        lines.append(f"    750nm : {imaging_params.laser_750_power}")
        lines.append("")
        lines.append("  Exposure Times (ms):")
        lines.append(f"    488nm : {imaging_params.exposure_488}")
        lines.append(f"    561nm : {imaging_params.exposure_561}")
        lines.append(f"    640nm : {imaging_params.exposure_640}")
        lines.append(f"    750nm : {imaging_params.exposure_750}")
        lines.append("")

    return "\n".join(lines)


def _generate_drift_plot(drift_df: pd.DataFrame, output_path: Path) -> bool:
    """Create a scatter plot of FOV positions coloured by round.

    Returns ``True`` if the plot was saved successfully, ``False`` if
    matplotlib is unavailable or the data is empty.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    if drift_df.empty:
        return False

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    rounds = sorted(drift_df["round"].unique())
    cmap = plt.cm.get_cmap("tab10", max(len(rounds), 1))

    # Left panel: delta_x vs delta_y scatter per round
    ax = axes[0]
    for idx, rnd in enumerate(rounds):
        subset = drift_df[drift_df["round"] == rnd]
        ax.scatter(
            subset["delta_x"],
            subset["delta_y"],
            c=[cmap(idx)],
            label=f"Round {rnd}",
            alpha=0.7,
            s=20,
            edgecolors="none",
        )
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.set_xlabel("delta_x")
    ax.set_ylabel("delta_y")
    ax.set_title("Per-FOV Drift (relative to reference round)")
    ax.legend(fontsize=7, loc="best")
    ax.set_aspect("equal", adjustable="datalim")

    # Right panel: displacement distribution per round (box plot)
    ax2 = axes[1]
    round_data = [
        drift_df[drift_df["round"] == rnd]["displacement"].values for rnd in rounds
    ]
    bp = ax2.boxplot(round_data, labels=[str(r) for r in rounds], patch_artist=True)
    for idx, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(cmap(idx))
        patch.set_alpha(0.7)
    ax2.set_xlabel("Round")
    ax2.set_ylabel("Displacement")
    ax2.set_title("Displacement Distribution per Round")

    fig.suptitle("Position Drift Analysis", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return True


# ---------------------------------------------------------------------------
# Stage implementation
# ---------------------------------------------------------------------------

@register_stage("inspect_positions")
class InspectPositionsStage(PipelineStage):
    """Analyse position drift across rounds and parse microscope logs."""

    description = "Position drift analysis and microscope-log inspection"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that the standardized positions file from the index stage exists."""
        errors: list[str] = []

        positions_path = self._positions_path()
        if not positions_path.exists():
            errors.append(
                f"Normalised position file not found: {positions_path}  "
                "(has the 'index' stage been run?)"
            )

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if the drift report and summary already exist."""
        out = self.get_output_dir()
        return (out / "drift_report.csv").exists() and (out / "drift_summary.txt").exists()

    def run(self, dry_run: bool = False) -> StageResult:
        """Execute the inspect_positions analysis.

        Steps
        -----
        1. Read normalised position data from the ``index`` stage output.
        2. Compute per-FOV drift relative to the first round.
        3. Optionally parse the microscope log for imaging parameters.
        4. Write ``drift_report.csv``, ``drift_summary.txt``, and
           (if matplotlib is available) ``drift_plot.png``.
        """
        start_time = datetime.now()
        output_dir = self.get_output_dir()

        if dry_run:
            self.logger.info("[DRY RUN] Would write outputs to %s", output_dir)
            return StageResult(status="skipped", metadata={"dry_run": True})

        output_dir.mkdir(parents=True, exist_ok=True)

        # ---- 1. Read positions ----
        positions_path = self._positions_path()
        self.logger.info("Reading positions from %s", positions_path)

        try:
            positions_df = read_sheet(positions_path)
        except Exception as exc:
            return StageResult(
                status="failed",
                error=f"Failed to read position file: {exc}",
            )

        if positions_df.empty:
            self.logger.warning("Position file is empty -- nothing to analyse.")
            return StageResult(
                status="completed",
                output_files=[],
                metadata={"warning": "empty position file"},
            )

        # ---- 2. Compute drift report ----
        rounds_to_check = self.config.inspect_positions.rounds_to_check
        self.logger.info(
            "Computing drift (rounds_to_check=%s) ...",
            rounds_to_check if rounds_to_check else "all",
        )

        try:
            drift_df = _compute_drift_report(positions_df, rounds_to_check)
        except ValueError as exc:
            return StageResult(status="failed", error=str(exc))

        report_path = output_dir / "drift_report.csv"
        drift_df.to_csv(report_path, index=False)
        self.logger.info("Wrote drift report: %s (%d rows)", report_path, len(drift_df))

        # ---- 3. Optionally parse microscope log ----
        imaging_params: Optional[_ImagingParams] = None
        log_path = self._resolve_log_path()

        if log_path is not None and log_path.exists():
            self.logger.info("Parsing microscope log: %s", log_path)
            imaging_params = _parse_microscope_log(log_path)
        elif log_path is not None:
            self.logger.warning("Microscope log not found: %s", log_path)
        else:
            self.logger.info("No microscope log configured -- skipping log parse.")

        # ---- 4. Write summary ----
        summary_text = _format_drift_summary(drift_df, imaging_params)
        summary_path = output_dir / "drift_summary.txt"
        summary_path.write_text(summary_text, encoding="utf-8")
        self.logger.info("Wrote drift summary: %s", summary_path)

        output_files = [str(report_path), str(summary_path)]

        # ---- 5. Optional plot ----
        plot_path = output_dir / "drift_plot.png"
        if _generate_drift_plot(drift_df, plot_path):
            output_files.append(str(plot_path))
            self.logger.info("Wrote drift plot: %s", plot_path)
        else:
            self.logger.info(
                "Drift plot skipped (matplotlib unavailable or no data)."
            )

        # ---- Metadata ----
        metadata: dict[str, Any] = {
            "n_rounds": int(drift_df["round"].nunique()) if not drift_df.empty else 0,
            "n_fovs": int(drift_df["fov"].nunique()) if not drift_df.empty else 0,
            "mean_displacement": float(drift_df["displacement"].mean()) if not drift_df.empty else 0.0,
            "max_displacement": float(drift_df["displacement"].max()) if not drift_df.empty else 0.0,
        }

        if imaging_params is not None:
            metadata["z_stacks"] = imaging_params.z_stacks
            metadata["upper_planes"] = imaging_params.upper_planes
            metadata["lower_planes"] = imaging_params.lower_planes

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata=metadata,
        )

        self.write_run_metadata(result, start_time)

        return result

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _positions_path(self) -> Path:
        """Return the expected path to the standardized positions file."""
        return Path(self.config.paths.output_dir) / "index" / "positions.standardized.csv"

    def _resolve_log_path(self) -> Optional[Path]:
        """Determine the microscope log path from config.

        Resolution order:
        1. ``config.inspect_positions.log_file`` (explicit override)
        2. ``<raw_data_dir> / config.microscope.log_file_name`` (convention)
        """
        explicit = self.config.inspect_positions.log_file
        if explicit is not None:
            return Path(explicit)

        # Fall back to default log filename in the raw data directory.
        try:
            default_name = self.config.microscope.log_file_name
        except AttributeError:
            return None

        if default_name:
            return Path(self.config.paths.raw_data_dir) / default_name

        return None
