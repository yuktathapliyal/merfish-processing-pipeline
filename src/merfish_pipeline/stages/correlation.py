"""``correlation`` stage -- barcode correlation analysis against bulk RNA-seq data.

For each distance threshold this stage:

1. Loads the codebook CSV (maps barcode_id to gene_symbol).
2. Loads the barcodes CSV produced by MERlin (or the ``filter_barcodes``
   stage) which contains a ``mean_distance`` column.
3. Loads a bulk RNA-seq expression CSV (gene_symbol to TPM values).
4. Filters barcodes by ``mean_distance < threshold``, counts occurrences per
   barcode_id, merges with the codebook and bulk expression data, and computes
   Pearson and Spearman correlations on log2-transformed values.
5. Generates scatter plots with gene labels and saves them as a multi-page PDF.
6. Writes an ``info_{xp}.csv`` summary table with metrics across thresholds.

Algorithm
---------
For every configured distance threshold (default
``[0.5167, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25]``):

    a. Filter barcodes where ``mean_distance < threshold``.
    b. Count occurrences per ``barcode_id``.
    c. Merge counts with the codebook to obtain ``gene_symbol``.
    d. Identify blank barcodes via regex ``^[Bb]lank[-_]?\\d+$``.
    e. Merge with bulk expression data on ``gene_symbol``.
    f. Compute Pearson and Spearman correlations (``scipy.stats``).
    g. Log-transform: ``np.log2(count + 1)`` / ``np.log2(tpm + 1)``.
    h. Scatter plot with gene labels saved to a multi-page PDF.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")  # non-interactive backend -- safe for headless servers
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

# ---------------------------------------------------------------------------
# Blank barcode detection
# ---------------------------------------------------------------------------

#: Regex matching blank / control barcode names such as ``Blank_01``,
#: ``blank-3``, ``Blank02``, etc.
_BLANK_RE = re.compile(r"^[Bb]lank[-_]?\d+$")


def _is_blank(gene_symbol: str) -> bool:
    """Return ``True`` when *gene_symbol* matches the blank barcode pattern."""
    return bool(_BLANK_RE.match(str(gene_symbol)))


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _filter_distance(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Filter barcodes where ``mean_distance < threshold``.

    Parameters
    ----------
    df:
        Barcodes dataframe with a ``mean_distance`` column.
    threshold:
        Maximum mean distance (exclusive) to retain.

    Returns
    -------
    pd.DataFrame
        Filtered copy of the input dataframe.
    """
    return df.loc[df["mean_distance"] < threshold].copy()


def _count_barcodes(df: pd.DataFrame) -> pd.DataFrame:
    """Group by ``barcode_id`` and count occurrences.

    Returns a dataframe with columns ``barcode_id`` and ``counts``.
    """
    return df.groupby("barcode_id").size().reset_index(name="counts")


def _merge_codebook(
    counts_df: pd.DataFrame,
    codebook_df: pd.DataFrame,
) -> pd.DataFrame:
    """Inner-join *counts_df* with the codebook on ``barcode_id``.

    The codebook is expected to have a ``gene_symbol`` column (or ``name``
    which is renamed) and either an explicit ``barcode_id`` column or the
    row index is used as barcode_id (matching the reference scripts).
    """
    cb = codebook_df.copy()

    # The reference scripts use the codebook index as barcode_id and rename
    # the ``name`` column to ``gene_symbol`` / ``gene_name``.
    if "barcode_id" not in cb.columns:
        cb["barcode_id"] = cb.index

    if "gene_symbol" not in cb.columns:
        if "name" in cb.columns:
            cb = cb.rename(columns={"name": "gene_symbol"})
        elif "gene_name" in cb.columns:
            cb = cb.rename(columns={"gene_name": "gene_symbol"})

    return pd.merge(cb, counts_df, how="inner", on="barcode_id")


def _merge_bulk(
    merged_df: pd.DataFrame,
    bulk_df: pd.DataFrame,
) -> pd.DataFrame:
    """Inner-join *merged_df* with bulk expression data on ``gene_symbol``.

    The bulk expression dataframe is expected to contain a ``gene_symbol``
    column and a numeric expression column (``bulk_exp``, ``FPKM``, or
    ``TPM``).  The expression column is normalized to ``bulk_exp`` in the
    returned dataframe.
    """
    bk = bulk_df.copy()

    # Normalize expression column name to ``bulk_exp``.
    for candidate in ("bulk_exp", "FPKM", "TPM", "tpm", "fpkm"):
        if candidate in bk.columns and candidate != "bulk_exp":
            bk = bk.rename(columns={candidate: "bulk_exp"})
            break

    # Normalize gene column name.
    if "gene_symbol" not in bk.columns:
        for candidate in ("gene_name", "Gene", "gene"):
            if candidate in bk.columns:
                bk = bk.rename(columns={candidate: "gene_symbol"})
                break

    return pd.merge(bk, merged_df, how="inner", on="gene_symbol")


def _compute_correlation(
    df: pd.DataFrame,
) -> tuple[float, float, float, float]:
    """Compute Pearson and Spearman correlations on log-transformed values.

    The dataframe must contain ``log_counts`` and ``log_tpm`` columns.

    Returns
    -------
    tuple
        ``(pearson_r, pearson_p, spearman_r, spearman_p)``
    """
    if len(df) < 3:
        return (float("nan"), float("nan"), float("nan"), float("nan"))

    pr, pp = pearsonr(df["log_tpm"], df["log_counts"])
    sr, sp = spearmanr(df["log_tpm"], df["log_counts"])
    return float(pr), float(pp), float(sr), float(sp)


def _plot_correlation(
    df: pd.DataFrame,
    threshold: float,
    ax: plt.Axes,
    xp_name: str = "",
) -> None:
    """Draw a scatter plot of log-transformed counts vs bulk expression.

    Gene labels are placed next to each point.  Pearson and Spearman
    statistics are displayed in the upper-left corner of the plot.

    Parameters
    ----------
    df:
        Dataframe containing ``log_tpm``, ``log_counts``, and ``gene_symbol``.
    threshold:
        Distance threshold used for filtering (shown in the title).
    ax:
        Matplotlib axes on which to draw.
    xp_name:
        Experiment name used in the plot title.
    """
    ax.scatter(df["log_tpm"], df["log_counts"], s=20, alpha=0.7)

    # Add gene labels.
    for _, row in df.iterrows():
        ax.text(
            row["log_tpm"] + 0.02,
            row["log_counts"],
            str(row["gene_symbol"]),
            fontsize=5,
            alpha=0.8,
        )

    # Compute and annotate statistics.
    if len(df) >= 3:
        pr, _ = pearsonr(df["log_tpm"], df["log_counts"])
        sr, _ = spearmanr(df["log_tpm"], df["log_counts"])
        ax.text(
            0.01,
            0.95,
            f"Pearson = {pr:.2f}\nSpearman = {sr:.2f}",
            transform=ax.transAxes,
            verticalalignment="top",
            fontsize=10,
        )

    title = f"Bulk Correlation (dist < {threshold})"
    if xp_name:
        title = f"{xp_name} -- {title}"
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("log2(TPM + 1)", fontsize=11)
    ax.set_ylabel("log2(counts + 1)", fontsize=11)


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("correlation")
class CorrelationStage(PipelineStage):
    """Barcode correlation analysis against bulk RNA-seq expression data."""

    description = "Correlate decoded barcodes with bulk RNA-seq at varying distance thresholds"

    # Expected output filenames
    _PLOTS_NAME = "correlation_plots.pdf"
    _MERGED_DIR = "merged_counts"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that all required input files are present."""
        errors: list[str] = []

        barcodes_path = self._resolve_barcodes_file()
        if not barcodes_path.exists():
            errors.append(f"Barcodes file not found: {barcodes_path}")

        codebook_path = self._resolve_codebook_file()
        if codebook_path is None:
            errors.append("No codebook template configured (merlin.codebook_template)")
        elif not codebook_path.exists():
            errors.append(f"Codebook file not found: {codebook_path}")

        bulk_path = self._resolve_bulk_file()
        if bulk_path is None:
            errors.append("No bulk expression file configured (correlation.bulk_file)")
        elif not bulk_path.exists():
            errors.append(f"Bulk expression file not found: {bulk_path}")

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if the info CSV and plots PDF already exist."""
        out = self.get_output_dir()
        xp_name = self.config.experiment.name
        info_path = out / f"info_{xp_name}.csv"
        plots_path = out / self._PLOTS_NAME
        return info_path.exists() and plots_path.exists()

    def run(self, dry_run: bool = False) -> StageResult:
        """Execute the correlation analysis."""
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        output_dir.mkdir(parents=True, exist_ok=True)
        merged_dir = output_dir / self._MERGED_DIR
        merged_dir.mkdir(parents=True, exist_ok=True)

        xp_name = self.config.experiment.name
        thresholds = sorted(
            self.config.correlation.distance_thresholds, reverse=True
        )

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would compute correlation for %d thresholds", len(thresholds)
            )
            return StageResult(
                status="skipped",
                metadata={"dry_run": True, "thresholds": thresholds},
            )

        # ----------------------------------------------------------
        # Step 1: Load input data
        # ----------------------------------------------------------
        barcodes_path = self._resolve_barcodes_file()
        codebook_path = self._resolve_codebook_file()
        bulk_path = self._resolve_bulk_file()

        if codebook_path is None or bulk_path is None:
            return StageResult(
                status="failed",
                error="Codebook or bulk expression file not configured.",
            )

        self.logger.info("Loading barcodes from %s", barcodes_path)
        barcodes_df = read_sheet(barcodes_path)

        self.logger.info("Loading codebook from %s", codebook_path)
        codebook_df = read_sheet(codebook_path)

        self.logger.info("Loading bulk expression from %s", bulk_path)
        bulk_df = read_sheet(bulk_path)

        if "mean_distance" not in barcodes_df.columns:
            return StageResult(
                status="failed",
                error=(
                    "Barcodes file does not contain a 'mean_distance' column. "
                    f"Columns found: {list(barcodes_df.columns)}"
                ),
            )

        # ----------------------------------------------------------
        # Step 2: Iterate over distance thresholds
        # ----------------------------------------------------------
        info_rows: list[dict[str, Any]] = []
        output_files: list[str] = []

        plots_path = output_dir / self._PLOTS_NAME
        with PdfPages(str(plots_path)) as pdf:
            for threshold in thresholds:
                self.logger.info(
                    "Processing distance threshold %.4f ...", threshold
                )

                # 2a. Filter by mean_distance
                filtered_df = _filter_distance(barcodes_df, threshold)
                n_barcodes = len(filtered_df)
                self.logger.info(
                    "  Barcodes after filtering: %d", n_barcodes
                )

                if n_barcodes == 0:
                    self.logger.warning(
                        "  No barcodes pass threshold %.4f, skipping.", threshold
                    )
                    info_rows.append(
                        {
                            "threshold": threshold,
                            "n_barcodes": 0,
                            "n_genes": 0,
                            "n_blanks": 0,
                            "pearson_r": float("nan"),
                            "pearson_p": float("nan"),
                            "spearman_r": float("nan"),
                            "spearman_p": float("nan"),
                        }
                    )
                    continue

                # 2b. Count per barcode_id
                counts_df = _count_barcodes(filtered_df)

                # 2c. Merge with codebook
                merged_cb = _merge_codebook(counts_df, codebook_df)
                n_genes = len(merged_cb)

                # 2d. Identify blank barcodes
                blank_mask = merged_cb["gene_symbol"].apply(_is_blank)
                n_blanks = int(blank_mask.sum())
                blank_counts = int(
                    merged_cb.loc[blank_mask, "counts"].sum()
                ) if n_blanks > 0 else 0
                total_counts = int(merged_cb["counts"].sum())

                self.logger.info(
                    "  Genes: %d, Blanks: %d (counts: %d / %d)",
                    n_genes, n_blanks, blank_counts, total_counts,
                )

                # 2e. Merge with bulk expression (non-blank only)
                non_blank_df = merged_cb.loc[~blank_mask].copy()
                merged_bulk = _merge_bulk(non_blank_df, bulk_df)

                if merged_bulk.empty:
                    self.logger.warning(
                        "  No genes matched between codebook and bulk expression "
                        "at threshold %.4f.",
                        threshold,
                    )
                    info_rows.append(
                        {
                            "threshold": threshold,
                            "n_barcodes": total_counts,
                            "n_genes": n_genes,
                            "n_blanks": blank_counts,
                            "pearson_r": float("nan"),
                            "pearson_p": float("nan"),
                            "spearman_r": float("nan"),
                            "spearman_p": float("nan"),
                        }
                    )
                    continue

                # 2g. Log-transform
                merged_bulk["log_counts"] = np.log2(merged_bulk["counts"] + 1)
                merged_bulk["log_tpm"] = np.log2(merged_bulk["bulk_exp"] + 1)

                # Save per-threshold merged CSV
                threshold_csv = merged_dir / f"{xp_name}_{threshold}.csv"
                merged_bulk.to_csv(threshold_csv, index=False)
                output_files.append(str(threshold_csv))

                # 2f. Compute correlation
                pr, pp, sr, sp = _compute_correlation(merged_bulk)
                self.logger.info(
                    "  Pearson=%.4f (p=%.2e), Spearman=%.4f (p=%.2e)",
                    pr, pp, sr, sp,
                )

                info_rows.append(
                    {
                        "threshold": threshold,
                        "n_barcodes": total_counts,
                        "n_genes": n_genes,
                        "n_blanks": blank_counts,
                        "pearson_r": pr,
                        "pearson_p": pp,
                        "spearman_r": sr,
                        "spearman_p": sp,
                    }
                )

                # 2h. Generate scatter plot
                fig, ax = plt.subplots(figsize=(9, 9))
                _plot_correlation(merged_bulk, threshold, ax, xp_name=xp_name)
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

        output_files.append(str(plots_path))
        self.logger.info("Saved correlation plots to %s", plots_path)

        # ----------------------------------------------------------
        # Step 5: Write info table
        # ----------------------------------------------------------
        info_df = pd.DataFrame(info_rows)
        info_path = output_dir / f"info_{xp_name}.csv"
        info_df.to_csv(info_path, index=False)
        output_files.append(str(info_path))
        self.logger.info("Wrote info table to %s", info_path)

        # ----------------------------------------------------------
        # Result
        # ----------------------------------------------------------
        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata={
                "xp_name": xp_name,
                "n_thresholds": len(thresholds),
                "thresholds": thresholds,
                "total_input_barcodes": len(barcodes_df),
                "info_rows": len(info_rows),
            },
        )

        self.write_run_metadata(
            result,
            start_time,
            parameters={
                "distance_thresholds": thresholds,
                "barcodes_file": str(barcodes_path),
                "codebook_file": str(codebook_path),
                "bulk_file": str(bulk_path),
            },
        )

        return result

    # ------------------------------------------------------------------
    # Path resolution helpers
    # ------------------------------------------------------------------

    def _resolve_barcodes_file(self) -> Path:
        """Resolve the barcodes input file.

        Priority:
        1. Explicit ``correlation.barcodes_file`` config override.
        2. Filtered barcodes from the ``filter_barcodes`` stage output.
        3. MERlin ``ExportBarcodes/barcodes.csv`` under merlin_analysis.
        """
        output_dir = Path(self.config.paths.output_dir)

        # 1. Explicit config override.
        explicit = self.config.correlation.barcodes_file
        if explicit is not None:
            return Path(explicit)

        # 2. filter_barcodes stage output.
        filter_output = output_dir / "filter_barcodes" / "barcodes_filtered.csv"
        if filter_output.exists():
            return filter_output

        # 3. MERlin output.
        merlin_data_name = Path(self.config.paths.merlin_data_dir).name
        return (
            output_dir / "merlin_analysis" / merlin_data_name
            / "ExportBarcodes" / "barcodes.csv"
        )

    def _resolve_codebook_file(self) -> Path | None:
        """Resolve the codebook CSV path from the MERlin config."""
        if self.config.merlin.codebook_template is not None:
            return Path(self.config.merlin.codebook_template)
        return None

    def _resolve_bulk_file(self) -> Path | None:
        """Resolve the bulk expression CSV path from the correlation config."""
        if self.config.correlation.bulk_file is not None:
            return Path(self.config.correlation.bulk_file)
        return None
