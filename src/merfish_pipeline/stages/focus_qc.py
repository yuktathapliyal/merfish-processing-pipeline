"""``focus_qc`` stage -- per-FOV best-focus detection via Laplacian-of-Gaussian variance.

This stage wraps the logic from ``FindInFocusBeads.py`` (focus-score
computation) and ``InFocusSummary.py`` (heatmap + text summary) into a
single pipeline stage.

Algorithm
---------
1. Discover all bead-channel images under ``{raw_data_dir}/{bead_channel_folder}/``.
2. Parse each filename to extract imaging-round (IR), FOV, and z-slice indices.
3. Group images by ``(IR, FOV)`` and, for every z-slice in each group, compute
   a focus metric: ``Var(Laplacian(GaussianBlur(image)))``.
4. Select the z-slice with the highest score as the best-focus plane.
5. Write three output artefacts:

   - ``best_focus_slices.csv`` -- one row per FOV with best-z per IR.
   - ``heatmap.png``           -- colour-coded heatmap of best-z values.
   - ``summary.txt``           -- human-readable statistics on z-variation.
"""

from __future__ import annotations

import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import cv2
import matplotlib
matplotlib.use("Agg")  # non-interactive backend -- safe for headless servers
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from merfish_pipeline.io.path_utils import find_files_matching
from merfish_pipeline.io.tiff_io import read_tiff
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

# ---------------------------------------------------------------------------
# Filename parsing
# ---------------------------------------------------------------------------

#: Regex for the canonical merFISH filename pattern:
#: ``merFISH_{IR:2d}_{FOV:3d}_{z:2d}.(tif|tiff)``
_FILENAME_RE = re.compile(
    r"merFISH_(\d{2})_(\d{3})_(\d{2})\.(?:tif|tiff)$", re.IGNORECASE
)


def _parse_filename(path: Path) -> tuple[str, str, str]:
    """Extract ``(ir, fov, z)`` as zero-padded strings from a merFISH filename.

    Raises
    ------
    ValueError
        If the filename does not match the expected pattern.
    """
    m = _FILENAME_RE.match(path.name)
    if m is None:
        raise ValueError(f"Filename does not match merFISH pattern: {path.name}")
    return m.group(1), m.group(2), m.group(3)


# ---------------------------------------------------------------------------
# Focus scoring
# ---------------------------------------------------------------------------


def _focus_score(image: np.ndarray, sigma: float, ksize: int) -> float:
    """Compute Laplacian-of-Gaussian variance as a focus metric.

    Parameters
    ----------
    image:
        2-D grayscale image (uint8 or uint16).
    sigma:
        Standard deviation of the Gaussian pre-filter.
    ksize:
        Kernel size for the Laplacian operator (must be odd: 1, 3, 5, ...).

    Returns
    -------
    float
        Variance of the Laplacian response -- higher means sharper focus.
    """
    # Ensure 8-bit for cv2 compatibility when source is 16-bit.
    if image.dtype != np.uint8:
        # Normalise to 0-255 without saturating.
        img_f = image.astype(np.float32)
        lo, hi = img_f.min(), img_f.max()
        if hi - lo > 0:
            img_f = (img_f - lo) / (hi - lo) * 255.0
        image = img_f.astype(np.uint8)

    blurred = cv2.GaussianBlur(image, (0, 0), sigmaX=sigma, sigmaY=sigma)
    lap = cv2.Laplacian(blurred, cv2.CV_32F, ksize=ksize)
    return float(lap.var())


def _best_z_for_stack(
    files: list[Path], sigma: float, ksize: int
) -> tuple[str, list[float]]:
    """Find the best-focus z-slice within a single (IR, FOV) stack.

    Parameters
    ----------
    files:
        Paths to the z-slice images **sorted by z-index**.
    sigma:
        Gaussian sigma forwarded to :func:`_focus_score`.
    ksize:
        Laplacian kernel size forwarded to :func:`_focus_score`.

    Returns
    -------
    best_z:
        The zero-padded z-string of the slice with the highest score.
    scores:
        All computed scores in the same order as *files*.
    """
    scores: list[float] = []
    z_labels: list[str] = []

    for path in files:
        img = read_tiff(path)
        # Collapse multi-channel images to single plane if needed.
        if img.ndim == 3:
            img = img[0]
        scores.append(_focus_score(img, sigma=sigma, ksize=ksize))
        _, _, z = _parse_filename(path)
        z_labels.append(z)

    best_idx = int(np.argmax(scores))
    return z_labels[best_idx], scores


# ---------------------------------------------------------------------------
# Summary generation  (ported from InFocusSummary.py)
# ---------------------------------------------------------------------------


def _summarize(df: pd.DataFrame) -> dict[str, Any]:
    """Compute per-FOV and per-IR z-variation statistics.

    Parameters
    ----------
    df:
        DataFrame with a ``FOV`` column and one column per IR (``IR01``,
        ``IR02``, ...).  Values are integer best-z indices.
    """
    data = df.set_index("FOV")
    data = data.dropna(axis=1, how="all")

    # Convert remaining columns to numeric (they may arrive as strings).
    data = data.apply(pd.to_numeric, errors="coerce")

    z_min = data.min(axis=1)
    z_max = data.max(axis=1)
    z_range = z_max - z_min

    constant = z_range[z_range == 0].index.tolist()
    one_step = z_range[z_range == 1].index.tolist()
    varied = z_range[z_range > 1].index.tolist()

    const_values = {}
    for fov in constant:
        first_valid = data.loc[fov].dropna()
        const_values[fov] = int(first_valid.iloc[0]) if len(first_valid) > 0 else 0

    deviants: dict[str, dict[str, Any]] = {}
    for ir_col in data.columns:
        col = data[ir_col].dropna()
        if col.empty:
            continue
        mode_val = col.mode().iloc[0]
        bad_fovs = col.index[col != mode_val].tolist()
        deviants[ir_col] = {
            "mode": int(mode_val),
            "count": len(bad_fovs),
            "fovs": bad_fovs,
        }

    return {
        "constant": constant,
        "const_values": const_values,
        "one_step": one_step,
        "varied": varied,
        "deviants": deviants,
    }


def _write_summary(stats: dict[str, Any], path: Path) -> None:
    """Write a human-readable summary text file."""
    total_fovs = len(stats["constant"]) + len(stats["one_step"]) + len(stats["varied"])

    lines: list[str] = []
    lines.append(f"Total FOVs analyzed: {total_fovs}")
    lines.append("")

    # Section 1 -- constant FOVs
    const = stats["constant"]
    lines.append(
        f"1) FOVs with identical z-slice across all IRs "
        f"(max-min = 0) [count = {len(const)}]:"
    )
    if const:
        for fov in const:
            lines.append(f"   - FOV {fov}: slice {stats['const_values'][fov]}")
    else:
        lines.append("   (none)")
    lines.append("")

    # Section 2 -- one-step variation
    one = stats["one_step"]
    lines.append(
        f"2) FOVs with exactly 1-unit variation across IRs "
        f"(max-min = 1) [count = {len(one)}]:"
    )
    if one:
        lines.append("   " + ", ".join(f"FOV {f}" for f in one))
    else:
        lines.append("   (none)")
    lines.append("")

    # Section 3 -- larger variation
    var = stats["varied"]
    lines.append(
        f"3) FOVs with larger variation (max-min > 1) [count = {len(var)}]:"
    )
    if var:
        lines.append("   " + ", ".join(f"FOV {f}" for f in var))
    else:
        lines.append("   (none)")
    lines.append("")

    # Section 4 -- per-IR deviations
    lines.append("4) Per-IR deviations from mode:")
    for ir_col, info in stats["deviants"].items():
        entry = f"   - {ir_col}: mode = {info['mode']}, {info['count']} deviating FOV(s)"
        if info["fovs"]:
            entry += " -> " + ", ".join(f"FOV {f}" for f in info["fovs"])
        lines.append(entry)
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def _plot_heatmap(df: pd.DataFrame, path: Path) -> None:
    """Generate and save a heatmap of best-focus z-slices.

    Parameters
    ----------
    df:
        DataFrame with ``FOV`` column and one column per IR.
    path:
        Destination path for the PNG file.
    """
    data = df.set_index("FOV").dropna(axis=1, how="all")
    data = data.apply(pd.to_numeric, errors="coerce")
    n_rows, n_cols = data.shape

    cell_width = 0.45
    cell_height = 0.1
    fig_w = max(cell_width * n_cols, 3.0)
    fig_h = max(cell_height * n_rows, 3.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        data,
        cmap="Blues",
        annot=False,
        square=False,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={
            "label": "In focus z-slice",
            "shrink": 0.65,
            "aspect": 20,
        },
        ax=ax,
    )
    cbar = ax.collections[0].colorbar
    cbar.set_label("In focus z-slice", fontweight="bold")

    ax.set_title(
        "In focus Z-slice Heatmap by FOV & IR", fontweight="bold", pad=20
    )
    ax.set_xlabel("Imaging Round (IR)", fontweight="bold")
    ax.set_ylabel("FOV", fontweight="bold")
    ax.tick_params(axis="x", rotation=90, labelsize=10, pad=6)
    ax.tick_params(axis="y", labelsize=8)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.09)
    fig.savefig(str(path), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("focus_qc")
class FocusQCStage(PipelineStage):
    """Per-FOV focus quality check using Laplacian-of-Gaussian variance."""

    description = "Detect best-focus z-slice per FOV from bead-channel images"

    # Expected output filenames
    _CSV_NAME = "best_focus_slices.csv"
    _HEATMAP_NAME = "heatmap.png"
    _SUMMARY_NAME = "summary.txt"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that the bead-channel image directory exists and contains images."""
        errors: list[str] = []

        bead_dir = self._bead_dir()
        if not bead_dir.exists():
            errors.append(f"Bead channel directory does not exist: {bead_dir}")
        elif not bead_dir.is_dir():
            errors.append(f"Bead channel path is not a directory: {bead_dir}")
        else:
            # Quick check: at least one TIFF present.
            tiffs = list(bead_dir.glob("merFISH_*.[Tt][Ii][Ff][Ff]"))
            if not tiffs:
                errors.append(
                    f"No merFISH TIFF images found in bead channel directory: {bead_dir}"
                )

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if all three output artefacts already exist."""
        out = self.get_output_dir()
        return all(
            (out / name).exists()
            for name in (self._CSV_NAME, self._HEATMAP_NAME, self._SUMMARY_NAME)
        )

    def run(self, dry_run: bool = False) -> StageResult:
        """Execute focus detection, summary, and heatmap generation."""
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        output_dir.mkdir(parents=True, exist_ok=True)

        sigma = self.config.focus_qc.sigma
        ksize = self.config.focus_qc.ksize

        if dry_run:
            self.logger.info("[DRY RUN] Would analyse bead images in %s", self._bead_dir())
            return StageResult(
                status="skipped",
                metadata={"dry_run": True},
            )

        # ----------------------------------------------------------
        # Step 1: Discover bead-channel images
        # ----------------------------------------------------------
        bead_dir = self._bead_dir()
        self.logger.info("Scanning bead images in %s ...", bead_dir)

        files = find_files_matching(bead_dir, "merFISH_*.[Tt][Ii][Ff][Ff]")
        if not files:
            return StageResult(
                status="failed",
                error=f"No merFISH TIFF images found in {bead_dir}",
            )

        self.logger.info("Found %d bead-channel images.", len(files))

        # ----------------------------------------------------------
        # Step 2: Group by (IR, FOV)
        # ----------------------------------------------------------
        groups: dict[tuple[str, str], list[Path]] = defaultdict(list)
        parse_errors: list[str] = []

        for path in files:
            try:
                ir, fov, _z = _parse_filename(path)
                groups[(ir, fov)].append(path)
            except ValueError as exc:
                parse_errors.append(str(exc))

        if parse_errors:
            self.logger.warning(
                "%d file(s) skipped due to filename parse errors.", len(parse_errors)
            )

        if not groups:
            return StageResult(
                status="failed",
                error="No valid bead images after filename parsing.",
            )

        self.logger.info(
            "Grouped into %d (IR, FOV) stacks.", len(groups)
        )

        # ----------------------------------------------------------
        # Step 3: Compute best-focus z per (IR, FOV)
        # ----------------------------------------------------------
        results: dict[str, dict[str, str]] = defaultdict(dict)
        all_scores: dict[tuple[str, str], list[float]] = {}
        n_processed = 0

        for (ir, fov), path_list in sorted(groups.items()):
            sorted_stack = sorted(
                path_list, key=lambda p: int(_parse_filename(p)[2])
            )
            best_z, scores = _best_z_for_stack(sorted_stack, sigma=sigma, ksize=ksize)
            results[fov][ir] = best_z
            all_scores[(ir, fov)] = scores
            n_processed += 1
            if n_processed % 50 == 0:
                self.logger.info(
                    "  Processed %d / %d stacks ...", n_processed, len(groups)
                )

        self.logger.info("Focus scoring complete for %d stacks.", n_processed)

        # ----------------------------------------------------------
        # Step 4: Build output DataFrame
        # ----------------------------------------------------------
        fov_list = sorted(results.keys(), key=lambda x: int(x))
        # Discover all IRs actually present (sorted numerically).
        all_irs = sorted({ir for ir, _fov in groups.keys()}, key=lambda x: int(x))
        ir_columns = [f"IR{ir}" for ir in all_irs]

        rows: list[dict[str, Any]] = []
        for fov in fov_list:
            row: dict[str, Any] = {"FOV": fov}
            for ir in all_irs:
                row[f"IR{ir}"] = results[fov].get(ir, "")
            rows.append(row)

        df = pd.DataFrame(rows, columns=["FOV"] + ir_columns)

        # ----------------------------------------------------------
        # Step 5: Write outputs
        # ----------------------------------------------------------
        csv_path = output_dir / self._CSV_NAME
        df.to_csv(csv_path, index=False)
        self.logger.info("Wrote best-focus CSV: %s", csv_path)

        heatmap_path = output_dir / self._HEATMAP_NAME
        _plot_heatmap(df, heatmap_path)
        self.logger.info("Wrote heatmap: %s", heatmap_path)

        summary_path = output_dir / self._SUMMARY_NAME
        stats = _summarize(df)
        _write_summary(stats, summary_path)
        self.logger.info("Wrote summary: %s", summary_path)

        # ----------------------------------------------------------
        # Result
        # ----------------------------------------------------------
        output_files = [str(csv_path), str(heatmap_path), str(summary_path)]

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata={
                "n_fovs": len(fov_list),
                "n_irs": len(all_irs),
                "n_stacks": n_processed,
                "n_images": len(files),
                "sigma": sigma,
                "ksize": ksize,
                "parse_errors": len(parse_errors),
                "constant_fovs": len(stats["constant"]),
                "one_step_fovs": len(stats["one_step"]),
                "varied_fovs": len(stats["varied"]),
            },
        )

        self.write_run_metadata(result, start_time, parameters={"sigma": sigma, "ksize": ksize})

        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _bead_dir(self) -> Path:
        """Return the path to the bead-channel image folder."""
        return Path(self.config.paths.raw_data_dir) / self.config.raw_data.bead_channel_folder
