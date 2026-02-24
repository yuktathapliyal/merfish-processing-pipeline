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
    Three-panel figure: displacement trend, strip plot per round, and
    FOV x round heatmap.  Only generated when *matplotlib* is importable.
trajectory_plot.html  *(optional)*
    Interactive 3D stage-trajectory plot per z-slice (requires plotly).
    Each subplot shows one line per imaging round connecting all FOV
    positions (x, y, z) in tile order.  Hover shows FOV number and
    coordinates.  Defaults to the first 3 z-slices; configurable via
    ``inspect_positions.trajectory_z_slices``.
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
    lines.append(f"  Rounds analysed       : {[int(r) for r in sorted(drift_df['round'].unique())]}")
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
    """Create a 3-panel drift analysis figure.

    Panels:
    1. **Displacement trend** — median + IQR band over rounds.
    2. **Strip plot** — jittered per-FOV displacement by round.
    3. **FOV x Round heatmap** — displacement for every (FOV, round) pair.

    Returns ``True`` if the plot was saved successfully, ``False`` if
    matplotlib is unavailable or the data is empty.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
    except ImportError:
        return False

    if drift_df.empty:
        return False

    rounds = sorted(drift_df["round"].unique())
    fovs = sorted(drift_df["fov"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # --- Panel 1: Displacement trend (median + IQR) ---
    ax1 = axes[0]
    medians, q25s, q75s = [], [], []
    for rnd in rounds:
        vals = drift_df[drift_df["round"] == rnd]["displacement"]
        medians.append(vals.median())
        q25s.append(vals.quantile(0.25))
        q75s.append(vals.quantile(0.75))

    ax1.fill_between(rounds, q25s, q75s, alpha=0.25, color="#1f77b4", label="IQR (25-75%)")
    ax1.plot(rounds, medians, "o-", color="#1f77b4", linewidth=2, markersize=6, label="Median")
    ax1.set_xlabel("Round")
    ax1.set_ylabel("Displacement")
    ax1.set_title("Drift Trend Across Rounds")
    ax1.set_xticks(rounds)
    ax1.legend(fontsize=8)
    ax1.grid(axis="y", alpha=0.3)

    # --- Panel 2: Strip plot (jittered individual FOVs) ---
    ax2 = axes[1]
    rng = np.random.default_rng(42)
    cmap = plt.cm.get_cmap("tab10", max(len(rounds), 1))

    for idx, rnd in enumerate(rounds):
        subset = drift_df[drift_df["round"] == rnd]
        jitter = rng.uniform(-0.25, 0.25, size=len(subset))
        ax2.scatter(
            rnd + jitter,
            subset["displacement"],
            c=[cmap(idx)],
            alpha=0.5,
            s=12,
            edgecolors="none",
            label=f"Round {int(rnd)}",
        )

    # Overlay median markers
    ax2.plot(rounds, medians, "k_", markersize=18, markeredgewidth=2.5, zorder=5,
             label="Median")
    ax2.set_xlabel("Round")
    ax2.set_ylabel("Displacement")
    ax2.set_title("Per-FOV Displacement by Round")
    ax2.set_xticks(rounds)
    ax2.grid(axis="y", alpha=0.3)
    ax2.legend(fontsize=7, loc="upper right", ncol=2)

    # --- Panel 3: FOV x Round heatmap ---
    ax3 = axes[2]

    # Build a 2D array (FOV x Round)
    heatmap_data = np.full((len(fovs), len(rounds)), np.nan)
    fov_to_idx = {f: i for i, f in enumerate(fovs)}
    rnd_to_idx = {r: i for i, r in enumerate(rounds)}
    for _, row in drift_df.iterrows():
        fi = fov_to_idx.get(row["fov"])
        ri = rnd_to_idx.get(row["round"])
        if fi is not None and ri is not None:
            heatmap_data[fi, ri] = row["displacement"]

    vmax = np.nanpercentile(heatmap_data, 95)
    im = ax3.imshow(
        heatmap_data,
        aspect="auto",
        cmap="YlOrRd",
        norm=Normalize(vmin=0, vmax=vmax),
        interpolation="nearest",
    )
    ax3.set_xticks(range(len(rounds)))
    ax3.set_xticklabels([str(r) for r in rounds])
    ax3.set_xlabel("Round")
    ax3.set_ylabel(f"FOV (n={len(fovs)})")

    # Show FOV labels on y-axis only if manageable
    if len(fovs) <= 30:
        ax3.set_yticks(range(len(fovs)))
        ax3.set_yticklabels([str(f) for f in fovs], fontsize=6)
    else:
        # Show a subset of ticks
        step = max(1, len(fovs) // 10)
        tick_positions = list(range(0, len(fovs), step))
        ax3.set_yticks(tick_positions)
        ax3.set_yticklabels([str(fovs[i]) for i in tick_positions], fontsize=7)

    ax3.set_title("Displacement by FOV & Round")
    fig.colorbar(im, ax=ax3, label="Displacement", shrink=0.8)

    fig.suptitle("Position Drift Analysis", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return True


def _generate_trajectory_plot(
    positions_df: pd.DataFrame,
    output_path: Path,
    z_slices: Optional[list[int]] = None,
) -> bool:
    """Create an interactive 3D trajectory plot (HTML) per z-slice.

    For each z-slice, one 3D subplot shows one line per imaging round.
    Each line connects the (stage_pos_x, stage_pos_y, z_position) of all FOVs
    in tile_number order.  Overlapping lines = stable stage; diverging = drift.

    Hover shows FOV number, coordinates, and round for each point.

    Parameters
    ----------
    positions_df:
        Standardized position data with columns ``round``, ``tile_number``,
        ``stage_pos_x``, ``stage_pos_y``, and ``z_position_0 .. z_position_N``.
    output_path:
        Where to save the HTML file (should end in ``.html``).
    z_slices:
        Which z-slice indices to plot (e.g., ``[0, 1, 2]``).
        Defaults to the first 3 available z-position columns.

    Returns ``True`` if saved successfully.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        return False

    if positions_df.empty:
        return False

    # Discover available z-position columns
    z_cols = sorted(
        [c for c in positions_df.columns if c.startswith("z_position_")],
        key=lambda c: int(c.split("_")[-1]),
    )

    if not z_cols:
        return False

    # Determine which z-slices to plot
    if z_slices is not None:
        selected = [f"z_position_{i}" for i in z_slices if f"z_position_{i}" in z_cols]
    else:
        selected = z_cols[:3]  # default: first 3

    if not selected:
        return False

    n_panels = len(selected)
    rounds = sorted(positions_df["round"].unique())

    # plotly tab10-equivalent colors
    _COLORS = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ]

    fig = make_subplots(
        rows=1,
        cols=n_panels,
        subplot_titles=[f"Z-slice {int(c.split('_')[-1])}" for c in selected],
        specs=[[{"type": "scatter3d"}] * n_panels],
        horizontal_spacing=0.02,
    )

    for panel_idx, z_col in enumerate(selected):
        col = panel_idx + 1
        z_index = int(z_col.split("_")[-1])

        for rnd_idx, rnd in enumerate(rounds):
            rnd_df = positions_df[positions_df["round"] == rnd].sort_values("tile_number")

            if z_col not in rnd_df.columns or rnd_df[z_col].isna().all():
                continue

            xs = rnd_df["stage_pos_x"].values
            ys = rnd_df["stage_pos_y"].values
            zs = rnd_df[z_col].values
            fovs = rnd_df["tile_number"].values

            color = _COLORS[rnd_idx % len(_COLORS)]

            hover_text = [
                f"FOV {int(f)}<br>x={x:.3f}, y={y:.3f}<br>z={z:.3f}"
                for f, x, y, z in zip(fovs, xs, ys, zs)
            ]

            # Only show legend for the first panel to avoid duplicates
            show_legend = panel_idx == 0

            fig.add_trace(
                go.Scatter3d(
                    x=xs, y=ys, z=zs,
                    mode="lines+markers",
                    name=f"Round {int(rnd)}",
                    legendgroup=f"Round {int(rnd)}",
                    showlegend=show_legend,
                    line=dict(color=color, width=2),
                    marker=dict(color=color, size=3, opacity=0.8),
                    hovertext=hover_text,
                    hoverinfo="text",
                ),
                row=1, col=col,
            )

        # Axis labels
        fig.update_scenes(
            dict(
                xaxis_title="stage_pos_x",
                yaxis_title="stage_pos_y",
                zaxis_title=z_col,
            ),
            row=1, col=col,
        )

    fig.update_layout(
        title=dict(
            text="Stage Trajectory per Z-slice (line per round, points = FOVs)",
            font=dict(size=16),
        ),
        height=700,
        width=550 * n_panels,
        legend=dict(font=dict(size=11)),
    )

    fig.write_html(str(output_path), include_plotlyjs="cdn")

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
        """Check that position data is available.

        Accepts either the ``index`` stage's ``positions.standardized.csv``
        or per-round ``stagePos_Round#N.csv`` files produced by ``ims_convert``
        (ANDOR workflow).
        """
        errors: list[str] = []

        positions_path = self._positions_path()
        if not positions_path.exists():
            errors.append(
                f"No position data found. Looked for:\n"
                f"  1. {Path(self.config.paths.output_dir) / 'index' / 'positions.standardized.csv'}\n"
                f"  2. stagePos_Round#N.csv files in {self.config.paths.merlin_data_dir}\n"
                f"Has the 'index' or 'ims_convert' stage been run?"
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

        # ---- 5. Optional drift plot ----
        plot_path = output_dir / "drift_plot.png"
        if _generate_drift_plot(drift_df, plot_path):
            output_files.append(str(plot_path))
            self.logger.info("Wrote drift plot: %s", plot_path)
        else:
            self.logger.info(
                "Drift plot skipped (matplotlib unavailable or no data)."
            )

        # ---- 6. Optional 3D trajectory plot (interactive HTML) ----
        trajectory_path = output_dir / "trajectory_plot.html"
        z_slices = self.config.inspect_positions.trajectory_z_slices
        if _generate_trajectory_plot(positions_df, trajectory_path, z_slices):
            output_files.append(str(trajectory_path))
            self.logger.info("Wrote trajectory plot: %s", trajectory_path)
        else:
            self.logger.info(
                "Trajectory plot skipped (no z-position data or plotly unavailable)."
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
        """Return the path to the best available positions file.

        Resolution order:

        1. ``{output_dir}/index/positions.standardized.csv`` -- produced by the
           ``index`` stage (all microscope types).
        2. Per-round ``stagePos_Round#N.csv`` files in ``merlin_data_dir`` --
           produced by ``ims_convert`` (ANDOR workflow).  When found, they are
           merged into a single ``{output_dir}/inspect_positions/positions.merged.csv``
           with an added ``round`` column and that path is returned.
        """
        # Primary: index stage output
        index_path = Path(self.config.paths.output_dir) / "index" / "positions.standardized.csv"
        if index_path.exists():
            return index_path

        # Fallback: merge per-round stagePos CSVs from ims_convert
        merged = self._merge_ims_convert_positions()
        if merged is not None:
            return merged

        # Nothing found — return the primary path so validate_inputs can
        # report the standard "index stage not run" message.
        return index_path

    def _merge_ims_convert_positions(self) -> Path | None:
        """Merge per-round stagePos CSVs from ``ims_convert`` into one file.

        Returns the merged file path, or *None* if no stagePos files exist.
        """
        import re

        merlin_dir = Path(self.config.paths.merlin_data_dir)
        if not merlin_dir.is_dir():
            return None

        stage_pos_re = re.compile(r"stagePos_Round#(\d+)\.csv$")
        csv_files: list[tuple[int, Path]] = []
        for p in sorted(merlin_dir.iterdir()):
            m = stage_pos_re.match(p.name)
            if m:
                csv_files.append((int(m.group(1)), p))

        if not csv_files:
            return None

        frames: list[pd.DataFrame] = []
        for round_num, csv_path in csv_files:
            try:
                df = pd.read_csv(csv_path)
                df["round"] = round_num
                frames.append(df)
            except Exception as exc:
                self.logger.warning(
                    "Could not read stagePos file %s: %s", csv_path, exc
                )

        if not frames:
            return None

        merged_df = pd.concat(frames, ignore_index=True)
        # Reorder so 'round' is the first column
        cols = ["round"] + [c for c in merged_df.columns if c != "round"]
        merged_df = merged_df[cols]

        out_dir = self.get_output_dir()
        out_dir.mkdir(parents=True, exist_ok=True)
        merged_path = out_dir / "positions.merged.csv"
        merged_df.to_csv(merged_path, index=False)
        self.logger.info(
            "Merged %d stagePos files from ims_convert → %s (%d rows)",
            len(csv_files), merged_path, len(merged_df),
        )
        return merged_path

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
