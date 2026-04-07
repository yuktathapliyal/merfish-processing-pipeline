"""``barcode_qc`` stage -- post-MERlin QC report with barcode metrics.

Computes summary statistics and generates diagnostic plots from the decoded
barcodes CSV, codebook, and (optionally) MERlin's PlotPerformance output.

Outputs
-------
- ``{output_dir}/barcode_qc/qc_summary.csv``
- ``{output_dir}/barcode_qc/per_fov_stats.csv``
- ``{output_dir}/barcode_qc/per_gene_stats.csv``
- ``{output_dir}/barcode_qc/per_bit_stats.csv``  (only if intensity columns present)
- ``{output_dir}/barcode_qc/per_cell_stats.csv``  (only if Cell_ID present)
- ``{output_dir}/barcode_qc/qc_report.pdf``
- ``{output_dir}/barcode_qc/spatial_plots/{experiment}_FOV_{NNN}.pdf``  (per-FOV)
- ``{output_dir}/barcode_qc/run_metadata.json``
"""

from __future__ import annotations

import logging
import re
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

from merfish_pipeline.io.codebook import is_blank, load_codebook
from merfish_pipeline.io.columns import (
    LOCAL_X_CANDIDATES,
    LOCAL_Y_CANDIDATES,
    Z_CANDIDATES,
    detect_column,
)
from merfish_pipeline.io.sheet_io import read_sheet, write_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)


def _find_intensity_columns(df: pd.DataFrame) -> list[str]:
    """Find ``intensity_0``, ``intensity_1``, ... columns in *df*."""
    cols = sorted(
        [c for c in df.columns if re.match(r"^intensity_\d+$", c)],
        key=lambda c: int(c.split("_")[1]),
    )
    return cols


def _get_bit_columns(codebook: pd.DataFrame) -> list[str]:
    """Identify binary 0/1 bit columns in a codebook (e.g. RS0015, RS0083)."""
    skip = {"name", "gene_name", "gene_symbol", "id", "barcode_id"}
    candidates = [c for c in codebook.columns if c not in skip]
    bit_cols = []
    for c in candidates:
        vals = codebook[c].dropna().unique()
        if len(vals) > 0 and set(vals).issubset({0, 1, 0.0, 1.0}):
            bit_cols.append(c)
    return bit_cols


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------


def _compute_gene_stats(
    barcodes: pd.DataFrame, codebook: pd.DataFrame
) -> pd.DataFrame:
    """Per-gene barcode counts with blank flag."""
    counts = barcodes.groupby("barcode_id").size().reset_index(name="count")
    merged = pd.merge(
        codebook[["barcode_id", "gene_symbol"]],
        counts,
        on="barcode_id",
        how="left",
    )
    merged["count"] = merged["count"].fillna(0).astype(int)
    merged["is_blank"] = merged["gene_symbol"].apply(is_blank)
    return merged.sort_values("count", ascending=False).reset_index(drop=True)


def _compute_fov_stats(barcodes: pd.DataFrame, fov_col: str) -> pd.DataFrame:
    """Per-FOV barcode counts and mean intensity."""
    groups = barcodes.groupby(fov_col)
    stats = groups.size().reset_index(name="n_barcodes")
    if "mean_intensity" in barcodes.columns:
        intensity = groups["mean_intensity"].mean().reset_index(name="mean_intensity")
        stats = pd.merge(stats, intensity, on=fov_col)
    return stats


def _compute_cell_stats(barcodes: pd.DataFrame) -> pd.DataFrame | None:
    """Per-cell stats. Returns None if Cell_ID not present or all null."""
    if "Cell_ID" not in barcodes.columns:
        return None
    assigned = barcodes[barcodes["Cell_ID"].notna()].copy()
    if assigned.empty:
        return None

    groups = assigned.groupby("Cell_ID")
    stats = groups.size().reset_index(name="n_barcodes")

    # Count unique genes per cell via codebook barcode_id
    if "barcode_id" in assigned.columns:
        genes = groups["barcode_id"].nunique().reset_index(name="n_genes")
        stats = pd.merge(stats, genes, on="Cell_ID")

    return stats


def _compute_summary(
    barcodes: pd.DataFrame,
    gene_stats: pd.DataFrame,
    fov_stats: pd.DataFrame,
    cell_stats: pd.DataFrame | None,
    perf_df: pd.DataFrame | None,
    fov_col: str,
) -> dict:
    """Compute single-row QC summary metrics."""
    s: dict = {}

    # Barcode totals
    s["total_barcodes"] = len(barcodes)
    s["unique_genes"] = int(gene_stats[~gene_stats["is_blank"]]["gene_symbol"].nunique())
    s["unique_fovs"] = int(barcodes[fov_col].nunique())

    # Per-gene stats
    coding = gene_stats[~gene_stats["is_blank"]]["count"]
    s["barcodes_per_gene_mean"] = round(coding.mean(), 2) if len(coding) else 0
    s["barcodes_per_gene_std"] = round(coding.std(), 2) if len(coding) else 0
    s["barcodes_per_gene_cv"] = (
        round(coding.std() / coding.mean(), 3) if coding.mean() > 0 else 0
    )

    # Blanks
    blank_count = int(gene_stats.loc[gene_stats["is_blank"], "count"].sum())
    s["blank_barcode_count"] = blank_count
    s["blank_barcode_pct"] = (
        round(blank_count / len(barcodes) * 100, 2) if len(barcodes) > 0 else 0
    )

    # Per-FOV
    s["barcodes_per_fov_mean"] = round(fov_stats["n_barcodes"].mean(), 2)
    s["barcodes_per_fov_std"] = round(fov_stats["n_barcodes"].std(), 2)
    s["barcodes_per_fov_cv"] = (
        round(fov_stats["n_barcodes"].std() / fov_stats["n_barcodes"].mean(), 3)
        if fov_stats["n_barcodes"].mean() > 0
        else 0
    )

    # Intensity
    if "mean_intensity" in barcodes.columns:
        s["mean_intensity_median"] = round(barcodes["mean_intensity"].median(), 4)
    if "mean_distance" in barcodes.columns:
        s["mean_distance_median"] = round(barcodes["mean_distance"].median(), 4)

    # Per-cell (optional)
    if cell_stats is not None and len(cell_stats) > 0:
        bc = cell_stats["n_barcodes"]
        s["total_cells"] = len(cell_stats)
        s["barcodes_per_cell_mean"] = round(bc.mean(), 2)
        s["barcodes_per_cell_median"] = round(bc.median(), 2)
        s["barcodes_per_cell_std"] = round(bc.std(), 2)
        s["barcodes_per_cell_p25"] = round(bc.quantile(0.25), 2)
        s["barcodes_per_cell_p75"] = round(bc.quantile(0.75), 2)
        s["barcodes_per_cell_p95"] = round(bc.quantile(0.95), 2)
        if "n_genes" in cell_stats.columns:
            s["genes_per_cell_mean"] = round(cell_stats["n_genes"].mean(), 2)
            s["genes_per_cell_median"] = round(cell_stats["n_genes"].median(), 2)

        total_assigned = int(barcodes["Cell_ID"].notna().sum())
        s["pct_barcodes_in_cells"] = round(total_assigned / len(barcodes) * 100, 2)

    # PlotPerformance (optional)
    if perf_df is not None and len(perf_df) > 0:
        best_idx = perf_df["Pearson correlation"].idxmax()
        best = perf_df.loc[best_idx]
        s["optimal_threshold"] = float(best["Distance Threshold"])
        s["best_pearson"] = round(float(best["Pearson correlation"]), 4)
        s["best_spearman"] = round(float(best["Spearman correlation"]), 4)
        s["barcodes_at_optimal"] = int(best.get(
            "# detected barcodes", best.get("# detected barcodes (including blanks)", 0)
        ))
        s["blanks_at_optimal"] = int(best.get(
            "# detected control barcodes", best.get("# detected control (blank) barcodes", 0)
        ))

    return s


# ---------------------------------------------------------------------------
# SNR / error-rate metrics
# ---------------------------------------------------------------------------


def _compute_snr_stats(
    barcodes: pd.DataFrame,
    codebook: pd.DataFrame,
    intensity_cols: list[str],
) -> tuple[np.ndarray, dict] | None:
    """Compute per-barcode ON/OFF contrast ratio.

    ``contrast = (mean(ON) - mean(OFF)) / (mean(ON) + mean(OFF))``

    Returns ``(contrast_array, summary_dict)`` or *None* if computation is
    not possible.  ``contrast_array`` has one value per barcode row (NaN for
    barcodes with no codebook match).
    """
    bit_cols = _get_bit_columns(codebook)
    n_bits = len(intensity_cols)

    if not bit_cols:
        logger.warning("No bit columns found in codebook; skipping SNR.")
        return None

    usable = min(len(bit_cols), n_bits)
    if len(bit_cols) != n_bits:
        logger.warning(
            "Bit count mismatch: codebook=%d, intensity=%d; using %d.",
            len(bit_cols), n_bits, usable,
        )
    bit_cols = bit_cols[:usable]
    intensity_cols = intensity_cols[:usable]

    # Build mask matrix indexed by barcode_id
    max_bid = int(codebook["barcode_id"].max())
    mask_matrix = np.zeros((max_bid + 1, usable), dtype=bool)
    for _, row in codebook.iterrows():
        bid = int(row["barcode_id"])
        mask_matrix[bid] = [bool(row[c]) for c in bit_cols]

    intensities = barcodes[intensity_cols].values.astype(float)
    bid_col = barcodes["barcode_id"].values.astype(int)
    valid = (bid_col >= 0) & (bid_col <= max_bid)

    contrast = np.full(len(barcodes), np.nan)
    if valid.sum() == 0:
        logger.warning("No valid barcode_id matches; skipping SNR.")
        return None

    masks = mask_matrix[bid_col[valid]]  # (n_valid, usable)
    ints = intensities[valid]

    on_count = masks.sum(axis=1).astype(float)
    off_count = (~masks).sum(axis=1).astype(float)

    on_sum = np.where(masks, ints, 0.0).sum(axis=1)
    off_sum = np.where(~masks, ints, 0.0).sum(axis=1)

    mean_on = np.divide(on_sum, on_count, where=on_count > 0,
                        out=np.zeros_like(on_sum))
    mean_off = np.divide(off_sum, off_count, where=off_count > 0,
                         out=np.zeros_like(off_sum))

    denom = mean_on + mean_off
    c = np.divide(mean_on - mean_off, denom, where=denom > 0,
                  out=np.zeros_like(denom))
    contrast[valid] = c

    valid_contrast = contrast[~np.isnan(contrast)]
    summary = {
        "snr_contrast_median": round(float(np.median(valid_contrast)), 4),
        "snr_contrast_mean": round(float(np.mean(valid_contrast)), 4),
        "snr_contrast_std": round(float(np.std(valid_contrast)), 4),
    }
    return contrast, summary


def _compute_misid_rate(gene_stats: pd.DataFrame) -> float | None:
    """Misidentification rate estimated from blank barcodes.

    ``misid_rate = (blank_count / n_blank_genes) / (coding_count / n_coding_genes)``
    """
    blanks = gene_stats[gene_stats["is_blank"]]
    coding = gene_stats[~gene_stats["is_blank"]]

    n_blank_genes = len(blanks)
    n_coding_genes = len(coding)
    if n_blank_genes == 0 or n_coding_genes == 0:
        return None

    blank_count = float(blanks["count"].sum())
    coding_count = float(coding["count"].sum())
    if coding_count == 0:
        return None

    return round((blank_count / n_blank_genes) / (coding_count / n_coding_genes), 6)


def _compute_per_bit_stats(
    barcodes: pd.DataFrame,
    codebook: pd.DataFrame,
    intensity_cols: list[str],
) -> pd.DataFrame | None:
    """Per-bit intensity stats across all barcodes.

    For each bit position, partitions barcodes into ON-group and OFF-group
    using the codebook mask, then computes median intensity and contrast.
    """
    bit_cols = _get_bit_columns(codebook)
    usable = min(len(bit_cols), len(intensity_cols))
    if usable == 0:
        return None

    bit_cols = bit_cols[:usable]
    intensity_cols = intensity_cols[:usable]

    max_bid = int(codebook["barcode_id"].max())
    mask_matrix = np.zeros((max_bid + 1, usable), dtype=bool)
    for _, row in codebook.iterrows():
        bid = int(row["barcode_id"])
        mask_matrix[bid] = [bool(row[c]) for c in bit_cols]

    intensities = barcodes[intensity_cols].values.astype(float)
    bid_col = barcodes["barcode_id"].values.astype(int)
    valid = (bid_col >= 0) & (bid_col <= max_bid)
    masks = mask_matrix[bid_col[valid]]
    ints = intensities[valid]

    rows = []
    for i in range(usable):
        on_mask = masks[:, i]
        on_vals = ints[on_mask, i]
        off_vals = ints[~on_mask, i]

        med_on = float(np.median(on_vals)) if len(on_vals) > 0 else 0.0
        med_off = float(np.median(off_vals)) if len(off_vals) > 0 else 0.0
        denom = med_on + med_off
        contrast = (med_on - med_off) / denom if denom > 0 else 0.0

        rows.append({
            "bit_index": i,
            "bit_name": bit_cols[i],
            "median_on": round(med_on, 6),
            "median_off": round(med_off, 6),
            "contrast": round(contrast, 4),
            "n_on": int(on_mask.sum()),
            "n_off": int((~on_mask).sum()),
        })
    return pd.DataFrame(rows)


def _compute_fov_misid_rate(
    barcodes: pd.DataFrame,
    codebook: pd.DataFrame,
    gene_stats: pd.DataFrame,
    fov_col: str,
) -> pd.DataFrame:
    """Per-FOV misidentification rate."""
    n_blank_genes = int(gene_stats["is_blank"].sum())
    n_coding_genes = int((~gene_stats["is_blank"]).sum())

    if n_blank_genes == 0 or n_coding_genes == 0:
        return pd.DataFrame(columns=[fov_col, "misid_rate"])

    blank_ids = set(
        codebook.loc[codebook["gene_symbol"].apply(is_blank), "barcode_id"]
    )
    # NOTE: name this ``blank_mask`` (not ``is_blank``) to avoid shadowing
    # the imported ``is_blank`` function from :mod:`merfish_pipeline.io.codebook`.
    blank_mask = barcodes["barcode_id"].isin(blank_ids)

    fov_blank = blank_mask.groupby(barcodes[fov_col]).sum().reset_index(
        name="_blank_count"
    )
    fov_total = barcodes.groupby(fov_col).size().reset_index(name="_total")
    fov = pd.merge(fov_blank, fov_total, on=fov_col)
    fov["_coding_count"] = fov["_total"] - fov["_blank_count"]
    fov["misid_rate"] = np.where(
        fov["_coding_count"] > 0,
        (fov["_blank_count"] / n_blank_genes)
        / (fov["_coding_count"] / n_coding_genes),
        np.nan,
    )
    fov["misid_rate"] = fov["misid_rate"].round(6)
    return fov[[fov_col, "misid_rate"]]


# ---------------------------------------------------------------------------
# PlotPerformance reader
# ---------------------------------------------------------------------------


def _read_plot_performance(perf_dir: Path) -> pd.DataFrame | None:
    """Read and merge MERlin PlotPerformance CSVs into a single DataFrame."""
    if not perf_dir.is_dir():
        return None

    csvs = sorted(perf_dir.glob("*.csv"))
    # Filter to threshold CSVs (exclude info_ files)
    csvs = [c for c in csvs if not c.name.startswith("info_")]
    if not csvs:
        return None

    frames = []
    for csv_path in csvs:
        try:
            df = read_sheet(csv_path)
            if "Distance Threshold" in df.columns:
                frames.append(df)
        except Exception:
            continue

    if not frames:
        return None

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.drop_duplicates(subset=["Distance Threshold"])
    combined = combined.sort_values("Distance Threshold").reset_index(drop=True)
    return combined


# ---------------------------------------------------------------------------
# Plot generation
# ---------------------------------------------------------------------------


def _generate_report(
    gene_stats: pd.DataFrame,
    fov_stats: pd.DataFrame,
    barcodes: pd.DataFrame,
    cell_stats: pd.DataFrame | None,
    perf_df: pd.DataFrame | None,
    top_n: int,
    output_path: Path,
    snr_contrast: np.ndarray | None = None,
    per_bit_stats: pd.DataFrame | None = None,
) -> None:
    """Generate a multi-panel QC report PDF."""
    has_cells = cell_stats is not None and len(cell_stats) > 0
    has_perf = perf_df is not None and len(perf_df) > 0
    has_snr = snr_contrast is not None
    has_bits = per_bit_stats is not None and len(per_bit_stats) > 0

    n_panels = (3 + int(has_perf) + int(has_cells) + 1  # +1 for top genes
                + int(has_snr) + int(has_bits))
    n_cols = 2
    n_rows = (n_panels + 1) // 2

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 4 * n_rows))
    axes = axes.flatten()
    panel = 0

    # Panel 1: Barcode abundance (sorted, log scale)
    ax = axes[panel]
    panel += 1
    coding = gene_stats[~gene_stats["is_blank"]].sort_values("count", ascending=False)
    blanks = gene_stats[gene_stats["is_blank"]].sort_values("count", ascending=False)
    ax.bar(range(len(coding)), coding["count"].values, color="steelblue", label="Coding")
    if len(blanks) > 0:
        ax.bar(
            range(len(coding), len(coding) + len(blanks)),
            blanks["count"].values,
            color="salmon",
            label="Blank",
        )
    ax.set_yscale("log")
    ax.set_xlabel("Gene rank")
    ax.set_ylabel("Barcode count")
    ax.set_title("Barcode abundance by gene")
    ax.legend(fontsize=8)

    # Panel 2: Intensity distribution
    ax = axes[panel]
    panel += 1
    if "mean_intensity" in barcodes.columns:
        vals = barcodes["mean_intensity"].dropna()
        vals_log = np.log10(vals[vals > 0])
        ax.hist(vals_log, bins=50, color="steelblue", edgecolor="none")
        ax.set_xlabel("log10(mean_intensity)")
        ax.set_ylabel("Count")
    ax.set_title("Intensity distribution")

    # Panel 3: Barcodes per FOV
    ax = axes[panel]
    panel += 1
    ax.bar(fov_stats.index, fov_stats["n_barcodes"].values, color="steelblue")
    ax.set_xlabel("FOV")
    ax.set_ylabel("Barcodes")
    ax.set_title("Barcodes per FOV")

    # Panel 4: Distance threshold curve (if available)
    if has_perf:
        ax = axes[panel]
        panel += 1
        ax.plot(
            perf_df["Distance Threshold"],
            perf_df["Pearson correlation"],
            "o-",
            label="Pearson",
            color="steelblue",
        )
        ax.plot(
            perf_df["Distance Threshold"],
            perf_df["Spearman correlation"],
            "s--",
            label="Spearman",
            color="darkorange",
        )
        best_idx = perf_df["Pearson correlation"].idxmax()
        best_thresh = perf_df.loc[best_idx, "Distance Threshold"]
        ax.axvline(best_thresh, color="red", linestyle=":", alpha=0.7, label=f"Best={best_thresh}")
        ax.set_xlabel("Distance threshold")
        ax.set_ylabel("Correlation")
        ax.set_title("Correlation vs distance threshold")
        ax.legend(fontsize=8)

    # Panel 5: Barcodes per cell (if available)
    if has_cells:
        ax = axes[panel]
        panel += 1
        ax.hist(
            cell_stats["n_barcodes"].values,
            bins=50,
            color="steelblue",
            edgecolor="none",
        )
        median_val = cell_stats["n_barcodes"].median()
        ax.axvline(median_val, color="red", linestyle="--", label=f"Median={median_val:.0f}")
        ax.set_xlabel("Barcodes per cell")
        ax.set_ylabel("Count")
        ax.set_title(f"Barcodes per cell (n={len(cell_stats)})")
        ax.legend(fontsize=8)

    # Panel 6: Top N genes
    ax = axes[panel]
    panel += 1
    top = coding.head(top_n)
    y_pos = np.arange(len(top))
    ax.barh(y_pos, top["count"].values, color="steelblue")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top["gene_symbol"].values, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Barcode count")
    ax.set_title(f"Top {top_n} genes")

    # Panel 7: SNR contrast distribution (if available)
    if has_snr:
        ax = axes[panel]
        panel += 1
        valid_snr = snr_contrast[~np.isnan(snr_contrast)]
        ax.hist(valid_snr, bins=60, color="steelblue", edgecolor="none")
        med = float(np.median(valid_snr))
        ax.axvline(med, color="red", linestyle="--", linewidth=1.5,
                    label=f"Median={med:.3f}")
        # Quality threshold lines
        for thresh, clr, lbl in [
            (0.3, "orange", "Poor<0.3"),
            (0.5, "gold", "Marginal<0.5"),
            (0.7, "green", "Good>0.7"),
        ]:
            ax.axvline(thresh, color=clr, linestyle=":", alpha=0.6, label=lbl)
        ax.set_xlabel("ON/OFF contrast ratio")
        ax.set_ylabel("Count")
        ax.set_title("SNR contrast distribution")
        ax.legend(fontsize=7, loc="upper left")

    # Panel 8: Per-bit contrast bar chart (if available)
    if has_bits:
        ax = axes[panel]
        panel += 1
        contrasts = per_bit_stats["contrast"].values
        colors = [
            "green" if c >= 0.7 else "gold" if c >= 0.5
            else "orange" if c >= 0.3 else "red"
            for c in contrasts
        ]
        ax.bar(range(len(contrasts)), contrasts, color=colors)
        ax.set_xticks(range(len(contrasts)))
        ax.set_xticklabels(per_bit_stats["bit_name"].values, rotation=45,
                           ha="right", fontsize=7)
        ax.set_ylabel("Contrast ratio")
        ax.set_title("Per-bit ON/OFF contrast")
        ax.set_ylim(0, 1.05)
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.4)

    # Hide unused axes
    for i in range(panel, len(axes)):
        axes[i].set_visible(False)

    fig.tight_layout()

    with PdfPages(output_path) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Per-FOV spatial scatter plots
# ---------------------------------------------------------------------------


def _generate_spatial_plots(
    barcodes: pd.DataFrame,
    fov_col: str,
    z_col: str,
    x_col: str,
    y_col: str,
    output_dir: Path,
    experiment_name: str,
    n_cols: int = 3,
) -> list[str]:
    """Generate per-FOV spatial scatter PDFs colored by ``mean_distance``.

    Each FOV gets its own PDF with a grid of subplots (one per Z-slice).
    Returns a list of output file paths.
    """
    if "mean_distance" not in barcodes.columns:
        logger.warning("Skipping spatial plots: 'mean_distance' column not found.")
        return []

    output_dir.mkdir(parents=True, exist_ok=True)
    output_files: list[str] = []

    # Global limits for consistent scales across FOVs
    vmax = float(barcodes["mean_distance"].max())
    xmax = float(barcodes[x_col].max())
    ymax = float(barcodes[y_col].max())
    offset = 25

    all_z = sorted(barcodes[z_col].unique())
    n_z = len(all_z)
    n_rows = -(-n_z // n_cols)  # ceil division

    for fov, df_fov in barcodes.groupby(fov_col):
        fig, axs = plt.subplots(
            n_rows, n_cols, figsize=(4.5 * n_cols, 4 * n_rows),
            constrained_layout=True,
        )
        axs = np.asarray(axs).flatten()

        for ax in axs:
            ax.set_visible(False)

        sc = None
        z_groups = df_fov.groupby(z_col)

        for z_idx, z_val in enumerate(all_z):
            if z_idx >= len(axs):
                break
            ax = axs[z_idx]
            ax.set_visible(True)

            if z_val in z_groups.groups:
                df_z = z_groups.get_group(z_val)
                sc = ax.scatter(
                    df_z[x_col].values,
                    df_z[y_col].values,
                    c=df_z["mean_distance"].values,
                    cmap="Reds_r",
                    vmin=0,
                    vmax=vmax,
                    s=3,
                    alpha=1,
                )
            ax.set_xlim(-offset, xmax + offset)
            ax.set_ylim(-offset, ymax + offset)
            ax.set_title(f"Z slice {int(z_val) + 1}")
            ax.set_aspect("equal", adjustable="datalim")

        if sc is not None:
            cbar = fig.colorbar(
                sc, ax=axs.tolist(), orientation="vertical",
            )
            cbar.set_label("Mean Distance to codebook")

        fig.suptitle(f"{experiment_name} -- FOV {int(fov):03d}", fontsize=14)

        pdf_path = output_dir / f"{experiment_name}_FOV_{int(fov):03d}.pdf"
        fig.savefig(pdf_path)
        plt.close(fig)
        output_files.append(str(pdf_path))

    return output_files


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("barcode_qc")
class BarcodeQCStage(PipelineStage):
    """Post-MERlin QC report with barcode metrics and diagnostic plots."""

    description = "Generate QC metrics and diagnostic plots from decoded barcodes"

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []

        barcodes_path = self._resolve_barcodes_path()
        if barcodes_path is None:
            errors.append(
                "Cannot locate barcodes CSV. Set barcode_qc.barcodes_file "
                "or ensure cell_assignment / filter_barcodes / MERlin output exists."
            )
        elif not barcodes_path.exists():
            errors.append(f"Barcodes file does not exist: {barcodes_path}")

        codebook_path = self._resolve_codebook_path()
        if codebook_path is None:
            errors.append("No codebook configured (merlin.codebook_template).")
        elif not codebook_path.exists():
            errors.append(f"Codebook file not found: {codebook_path}")

        return errors

    def check_outputs_exist(self) -> bool:
        out = self.get_output_dir()
        return (out / "qc_summary.csv").exists() and (out / "qc_report.pdf").exists()

    def run(self, dry_run: bool = False) -> StageResult:
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        cfg = self.config.barcode_qc

        barcodes_path = self._resolve_barcodes_path()
        codebook_path = self._resolve_codebook_path()

        if barcodes_path is None or codebook_path is None:
            return StageResult(
                status="failed",
                error="Could not resolve barcodes or codebook path.",
            )

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would generate QC report from %s with codebook %s",
                barcodes_path,
                codebook_path,
            )
            return StageResult(status="skipped", metadata={"dry_run": True})

        # Load data
        self.logger.info("Loading barcodes from %s ...", barcodes_path)
        barcodes = read_sheet(barcodes_path)
        self.logger.info("Loaded %d barcodes.", len(barcodes))

        codebook = load_codebook(codebook_path)

        # Detect FOV column
        fov_col = "fov"
        for c in ["fov", "FOV", "Fov"]:
            if c in barcodes.columns:
                fov_col = c
                break

        # Compute stats
        gene_stats = _compute_gene_stats(barcodes, codebook)
        fov_stats = _compute_fov_stats(barcodes, fov_col)
        cell_stats = _compute_cell_stats(barcodes)

        # Try to read PlotPerformance
        perf_df = None
        perf_dir = self._resolve_plot_performance_dir()
        if perf_dir is not None:
            perf_df = _read_plot_performance(perf_dir)
            if perf_df is not None:
                self.logger.info(
                    "Read PlotPerformance data: %d thresholds", len(perf_df)
                )

        summary = _compute_summary(
            barcodes, gene_stats, fov_stats, cell_stats, perf_df, fov_col
        )

        # --- SNR & error-rate metrics (optional) ---
        snr_contrast: np.ndarray | None = None
        per_bit_stats: pd.DataFrame | None = None
        intensity_cols = _find_intensity_columns(barcodes)

        if intensity_cols:
            snr_result = _compute_snr_stats(barcodes, codebook, intensity_cols)
            if snr_result is not None:
                snr_contrast, snr_summary = snr_result
                summary.update(snr_summary)
                self.logger.info(
                    "SNR contrast: median=%.4f, mean=%.4f",
                    snr_summary["snr_contrast_median"],
                    snr_summary["snr_contrast_mean"],
                )
                # Per-FOV SNR median
                barcodes["_snr_contrast"] = snr_contrast
                fov_snr = (
                    barcodes.dropna(subset=["_snr_contrast"])
                    .groupby(fov_col)["_snr_contrast"]
                    .median()
                    .reset_index(name="snr_contrast_median")
                )
                fov_stats = pd.merge(fov_stats, fov_snr, on=fov_col, how="left")
                barcodes.drop(columns=["_snr_contrast"], inplace=True)

            per_bit_stats = _compute_per_bit_stats(
                barcodes, codebook, intensity_cols
            )
        else:
            self.logger.info(
                "No per-bit intensity columns found; skipping SNR metrics."
            )

        # Misidentification rate (uses gene_stats, not intensity columns)
        misid = _compute_misid_rate(gene_stats)
        if misid is not None:
            summary["misid_rate"] = misid
            self.logger.info("Misidentification rate: %.4f%%", misid * 100)
            # Per-FOV misid_rate
            fov_mr = _compute_fov_misid_rate(
                barcodes, codebook, gene_stats, fov_col
            )
            if len(fov_mr) > 0:
                fov_stats = pd.merge(fov_stats, fov_mr, on=fov_col, how="left")

        # Write outputs
        output_dir.mkdir(parents=True, exist_ok=True)
        output_files: list[str] = []

        summary_path = output_dir / "qc_summary.csv"
        pd.DataFrame([summary]).to_csv(summary_path, index=False)
        output_files.append(str(summary_path))

        fov_path = output_dir / "per_fov_stats.csv"
        write_sheet(fov_stats, fov_path)
        output_files.append(str(fov_path))

        gene_path = output_dir / "per_gene_stats.csv"
        write_sheet(gene_stats, gene_path)
        output_files.append(str(gene_path))

        if cell_stats is not None:
            cell_path = output_dir / "per_cell_stats.csv"
            write_sheet(cell_stats, cell_path)
            output_files.append(str(cell_path))

        if per_bit_stats is not None:
            bit_path = output_dir / "per_bit_stats.csv"
            write_sheet(per_bit_stats, bit_path)
            output_files.append(str(bit_path))
            self.logger.info("Wrote per-bit stats: %s", bit_path)

        # Generate PDF report
        report_path = output_dir / "qc_report.pdf"
        _generate_report(
            gene_stats=gene_stats,
            fov_stats=fov_stats,
            barcodes=barcodes,
            cell_stats=cell_stats,
            perf_df=perf_df,
            top_n=cfg.top_n_genes,
            output_path=report_path,
            snr_contrast=snr_contrast,
            per_bit_stats=per_bit_stats,
        )
        output_files.append(str(report_path))
        self.logger.info("Wrote QC report: %s", report_path)

        # Generate per-FOV spatial scatter plots
        if cfg.spatial_plots_enabled:
            try:
                x_col = detect_column(barcodes, LOCAL_X_CANDIDATES, "x")
                y_col = detect_column(barcodes, LOCAL_Y_CANDIDATES, "y")
                z_col_detected = detect_column(barcodes, Z_CANDIDATES, "z")

                spatial_dir = output_dir / "spatial_plots"
                spatial_files = _generate_spatial_plots(
                    barcodes=barcodes,
                    fov_col=fov_col,
                    z_col=z_col_detected,
                    x_col=x_col,
                    y_col=y_col,
                    output_dir=spatial_dir,
                    experiment_name=self.config.experiment.name,
                    n_cols=cfg.spatial_plots_columns,
                )
                output_files.extend(spatial_files)
                self.logger.info(
                    "Generated %d spatial scatter plots in %s",
                    len(spatial_files), spatial_dir,
                )
            except ValueError as exc:
                self.logger.warning("Skipping spatial plots: %s", exc)

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata=summary,
        )

        self.write_run_metadata(result, start_time, parameters=summary)
        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_barcodes_path(self) -> Path | None:
        """Resolve barcodes CSV: explicit → cell_assignment → filter_barcodes → MERlin."""
        explicit = self.config.barcode_qc.barcodes_file
        if explicit is not None:
            return Path(explicit)

        output_dir = Path(self.config.paths.output_dir)

        # cell_assignment output (has Cell_ID)
        ca_output = output_dir / "cell_assignment" / "barcodes_assigned.csv"
        if ca_output.exists():
            return ca_output

        # filter_barcodes output
        fb_output = output_dir / "filter_barcodes" / "barcodes_filtered.csv"
        if fb_output.exists():
            return fb_output

        # MERlin output
        merlin_data_name = Path(self.config.paths.merlin_data_dir).name
        merlin_output = (
            output_dir / "merlin_analysis" / merlin_data_name
            / "ExportBarcodes" / "barcodes.csv"
        )
        if merlin_output.exists():
            return merlin_output

        return None

    def _resolve_codebook_path(self) -> Path | None:
        if self.config.merlin.codebook_template is not None:
            return Path(self.config.merlin.codebook_template)
        return None

    def _resolve_plot_performance_dir(self) -> Path | None:
        """Try to find MERlin's PlotPerformance directory."""
        output_dir = Path(self.config.paths.output_dir)
        merlin_data_name = Path(self.config.paths.merlin_data_dir).name
        candidate = output_dir / "merlin_analysis" / merlin_data_name / "PlotPerformance"
        if candidate.is_dir():
            return candidate
        return None
