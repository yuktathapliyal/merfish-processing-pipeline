"""``barcode_qc`` stage -- post-MERlin QC report with barcode metrics.

Computes summary statistics and generates diagnostic plots from the decoded
barcodes CSV, codebook, and (optionally) MERlin's PlotPerformance output.

Outputs
-------
- ``{output_dir}/barcode_qc/qc_summary.csv``
- ``{output_dir}/barcode_qc/per_fov_stats.csv``
- ``{output_dir}/barcode_qc/per_gene_stats.csv``
- ``{output_dir}/barcode_qc/per_cell_stats.csv``  (only if Cell_ID present)
- ``{output_dir}/barcode_qc/qc_report.pdf``
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

from merfish_pipeline.io.sheet_io import read_sheet, write_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)

_BLANK_RE = re.compile(r"^[Bb]lank[-_]?\d+$")


def _is_blank(gene_symbol: str) -> bool:
    return bool(_BLANK_RE.match(str(gene_symbol)))


# ---------------------------------------------------------------------------
# Codebook helpers (mirrors correlation.py pattern)
# ---------------------------------------------------------------------------


def _load_codebook(codebook_path: Path) -> pd.DataFrame:
    """Load codebook and normalise column names."""
    cb = read_sheet(codebook_path)
    if "barcode_id" not in cb.columns:
        cb["barcode_id"] = cb.index
    if "gene_symbol" not in cb.columns:
        if "name" in cb.columns:
            cb = cb.rename(columns={"name": "gene_symbol"})
        elif "gene_name" in cb.columns:
            cb = cb.rename(columns={"gene_name": "gene_symbol"})
    return cb


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
    merged["is_blank"] = merged["gene_symbol"].apply(_is_blank)
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
) -> None:
    """Generate a multi-panel QC report PDF."""
    has_cells = cell_stats is not None and len(cell_stats) > 0
    has_perf = perf_df is not None and len(perf_df) > 0

    n_panels = 3 + int(has_perf) + int(has_cells) + 1  # +1 for top genes
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

    # Hide unused axes
    for i in range(panel, len(axes)):
        axes[i].set_visible(False)

    fig.tight_layout()

    with PdfPages(output_path) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


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

        codebook = _load_codebook(codebook_path)

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
        )
        output_files.append(str(report_path))
        self.logger.info("Wrote QC report: %s", report_path)

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
