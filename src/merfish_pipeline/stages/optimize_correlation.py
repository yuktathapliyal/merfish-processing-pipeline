"""``optimize_correlation`` stage -- find gene subgroups maximizing bulk correlation.

Uses simulated annealing to find, for each target group size, the subset of
genes that maximizes Pearson correlation between merFISH barcode log-counts
and bulk RNA-seq log-TPM.

Reads the merged counts CSV produced by the ``correlation`` stage.

Outputs
-------
- ``{output_dir}/optimize_correlation/correlation_trend.csv``
- ``{output_dir}/optimize_correlation/correlation_trend.png``
- ``{output_dir}/optimize_correlation/optimal_genes.csv``
- ``{output_dir}/optimize_correlation/detailed_results.xlsx``
- ``{output_dir}/optimize_correlation/run_metadata.json``
"""

from __future__ import annotations

import logging
import math
import random
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Simulated annealing helpers
# ---------------------------------------------------------------------------


def _calculate_correlation(
    log_tpm: np.ndarray,
    log_counts: np.ndarray,
    indices: set[int],
) -> float:
    """Pearson r on the subset defined by *indices*. Returns -inf on failure."""
    idx = list(indices)
    if len(idx) < 2:
        return float("-inf")
    try:
        r, _ = pearsonr(log_tpm[idx], log_counts[idx])
        if np.isnan(r):
            return float("-inf")
        return float(r)
    except Exception:
        return float("-inf")


def _optimize_group(
    log_tpm: np.ndarray,
    log_counts: np.ndarray,
    n_genes: int,
    target_size: int,
    *,
    max_iterations: int,
    initial_temperature: float,
    cooling_rate: float,
    rng: random.Random,
) -> tuple[set[int], float]:
    """Single simulated-annealing run for a given target group size.

    Returns ``(best_indices, best_correlation)``.
    """
    # Initialize random group
    all_indices = list(range(n_genes))
    current = set(rng.sample(all_indices, min(target_size, n_genes)))
    current_r = _calculate_correlation(log_tpm, log_counts, current)

    best = set(current)
    best_r = current_r
    temperature = initial_temperature

    for _ in range(max_iterations):
        # Propose: swap one gene out, one gene in
        candidate = set(current)
        outside = [i for i in all_indices if i not in candidate]
        if not outside:
            break

        remove_gene = rng.choice(list(candidate))
        add_gene = rng.choice(outside)
        candidate.discard(remove_gene)
        candidate.add(add_gene)

        candidate_r = _calculate_correlation(log_tpm, log_counts, candidate)
        delta = candidate_r - current_r

        # Metropolis criterion
        if delta > 0:
            accept = True
        elif temperature > 1e-10:
            exponent = max(delta / temperature, -700)
            accept = rng.random() < math.exp(exponent)
        else:
            accept = False

        if accept:
            current = candidate
            current_r = candidate_r

        if current_r > best_r:
            best = set(current)
            best_r = current_r

        temperature *= cooling_rate

    return best, best_r


def _find_progressive_groups(
    df: pd.DataFrame,
    *,
    size_range: range,
    correlation_threshold: float,
    n_attempts: int,
    max_iterations: int,
    initial_temperature: float,
    cooling_rate: float,
    rng: random.Random,
) -> tuple[dict[int, tuple[set[int], float]], list[tuple[int, float]]]:
    """Run SA optimization across a range of target sizes.

    Returns
    -------
    results_by_size : dict
        ``{size: (best_indices, best_r)}`` for sizes where threshold was met.
    correlation_trend : list
        ``[(size, best_r), ...]`` for all tested sizes.
    """
    log_tpm = df["log_tpm"].values
    log_counts = df["log_counts"].values
    n_genes = len(df)

    results_by_size: dict[int, tuple[set[int], float]] = {}
    correlation_trend: list[tuple[int, float]] = []

    for target_size in size_range:
        if target_size > n_genes:
            logger.info(
                "  Skipping size %d (exceeds %d available genes)",
                target_size, n_genes,
            )
            break

        best_indices: set[int] | None = None
        best_r = float("-inf")

        for attempt in range(n_attempts):
            indices, r = _optimize_group(
                log_tpm, log_counts, n_genes, target_size,
                max_iterations=max_iterations,
                initial_temperature=initial_temperature,
                cooling_rate=cooling_rate,
                rng=rng,
            )
            if r > best_r:
                best_r = r
                best_indices = indices

        correlation_trend.append((target_size, best_r))
        logger.info("  Size %3d: best r = %.4f", target_size, best_r)

        if best_r >= correlation_threshold and best_indices is not None:
            results_by_size[target_size] = (best_indices, best_r)

    return results_by_size, correlation_trend


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def _plot_correlation_trend(
    trend_df: pd.DataFrame,
    threshold: float,
    output_path: Path,
) -> None:
    """Line plot of group size vs best Pearson correlation."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(
        trend_df["group_size"], trend_df["best_correlation"],
        "o-", color="steelblue", linewidth=2, markersize=5,
    )
    ax.axhline(
        threshold, color="red", linestyle="--", alpha=0.7,
        label=f"Threshold = {threshold}",
    )
    ax.set_xlabel("Gene Group Size", fontsize=12)
    ax.set_ylabel("Best Pearson Correlation", fontsize=12)
    ax.set_title("Correlation vs Gene Group Size (Simulated Annealing)", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_detailed_results(
    results_by_size: dict[int, tuple[set[int], float]],
    df: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write an Excel workbook with a summary sheet and per-size gene sheets."""
    summary_rows = []
    sheets: dict[str, pd.DataFrame] = {}

    for size in sorted(results_by_size):
        indices, r = results_by_size[size]
        gene_data = df.iloc[list(indices)][
            ["gene_symbol", "log_tpm", "log_counts"]
        ].copy()
        gene_data = gene_data.sort_values("gene_symbol").reset_index(drop=True)
        gene_data["correlation"] = r

        summary_rows.append({
            "group_size": size,
            "correlation": round(r, 4),
            "n_genes": len(indices),
        })

        sheet_name = f"Group_{size}"
        if len(sheet_name) > 31:
            sheet_name = sheet_name[:31]
        sheets[sheet_name] = gene_data

    if not summary_rows:
        # Write an empty workbook with just the summary
        summary_rows.append({
            "group_size": 0,
            "correlation": float("nan"),
            "n_genes": 0,
        })

    summary_df = pd.DataFrame(summary_rows)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        for name, sheet_df in sheets.items():
            sheet_df.to_excel(writer, sheet_name=name, index=False)


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("optimize_correlation")
class OptimizeCorrelationStage(PipelineStage):
    """Find gene subgroups maximizing bulk RNA-seq correlation via simulated annealing."""

    description = (
        "Simulated annealing optimization of gene subgroups "
        "for maximum correlation with bulk RNA-seq"
    )

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []

        merged_path = self._resolve_merged_counts_file()
        if merged_path is None:
            errors.append(
                "Cannot locate merged counts from correlation stage. "
                "Run the 'correlation' stage first, or set "
                "optimize_correlation.distance_threshold explicitly."
            )
        elif not merged_path.exists():
            errors.append(f"Merged counts file not found: {merged_path}")

        cfg = self.config.optimize_correlation
        if cfg.size_range_start >= cfg.size_range_end:
            errors.append(
                f"size_range_start ({cfg.size_range_start}) must be < "
                f"size_range_end ({cfg.size_range_end})"
            )
        if cfg.size_range_step < 1:
            errors.append("size_range_step must be >= 1")

        return errors

    def check_outputs_exist(self) -> bool:
        out = self.get_output_dir()
        return (
            (out / "correlation_trend.csv").exists()
            and (out / "optimal_genes.csv").exists()
        )

    def run(self, dry_run: bool = False) -> StageResult:
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        cfg = self.config.optimize_correlation

        merged_path = self._resolve_merged_counts_file()
        if merged_path is None or not merged_path.exists():
            return StageResult(
                status="failed",
                error="Could not resolve merged counts file.",
            )

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would optimize correlation from %s", merged_path,
            )
            return StageResult(status="skipped", metadata={"dry_run": True})

        # Load merged counts
        self.logger.info("Loading merged counts from %s", merged_path)
        df = read_sheet(merged_path)

        for required_col in ("log_tpm", "log_counts", "gene_symbol"):
            if required_col not in df.columns:
                return StageResult(
                    status="failed",
                    error=(
                        f"Merged counts missing required column '{required_col}'. "
                        f"Available: {list(df.columns)}"
                    ),
                )

        # Drop rows with NaN in critical columns
        df = df.dropna(subset=["log_tpm", "log_counts", "gene_symbol"]).reset_index(drop=True)
        n_genes = len(df)
        self.logger.info("Loaded %d genes for optimization.", n_genes)

        # Build size range, capped at available genes
        effective_end = min(cfg.size_range_end, n_genes) + 1
        size_range = range(cfg.size_range_start, effective_end, cfg.size_range_step)
        if not list(size_range):
            return StageResult(
                status="failed",
                error=(
                    f"Empty size range: start={cfg.size_range_start}, "
                    f"end={effective_end - 1}, step={cfg.size_range_step}, "
                    f"available genes={n_genes}"
                ),
            )

        # Initialize RNG
        rng = random.Random(cfg.random_seed)

        # Run progressive optimization
        self.logger.info(
            "Running SA optimization: sizes %d-%d (step %d), "
            "%d attempts, %d iterations each",
            cfg.size_range_start, effective_end - 1, cfg.size_range_step,
            cfg.n_attempts, cfg.max_iterations,
        )

        results_by_size, correlation_trend = _find_progressive_groups(
            df,
            size_range=size_range,
            correlation_threshold=cfg.correlation_threshold,
            n_attempts=cfg.n_attempts,
            max_iterations=cfg.max_iterations,
            initial_temperature=cfg.initial_temperature,
            cooling_rate=cfg.cooling_rate,
            rng=rng,
        )

        # Write outputs
        output_dir.mkdir(parents=True, exist_ok=True)
        output_files: list[str] = []

        # 1. Correlation trend CSV
        trend_df = pd.DataFrame(correlation_trend, columns=["group_size", "best_correlation"])
        trend_path = output_dir / "correlation_trend.csv"
        trend_df.to_csv(trend_path, index=False)
        output_files.append(str(trend_path))

        # 2. Correlation trend plot
        trend_plot_path = output_dir / "correlation_trend.png"
        _plot_correlation_trend(trend_df, cfg.correlation_threshold, trend_plot_path)
        output_files.append(str(trend_plot_path))

        # 3. Optimal genes CSV (best correlation group)
        optimal_path = output_dir / "optimal_genes.csv"
        if results_by_size:
            # Pick group with highest correlation
            best_size = max(results_by_size, key=lambda s: results_by_size[s][1])
            best_indices, best_r = results_by_size[best_size]
            optimal_df = df.iloc[list(best_indices)][
                ["gene_symbol", "log_tpm", "log_counts"]
            ].copy()
            optimal_df["correlation"] = best_r
            optimal_df["group_size"] = best_size
            optimal_df = optimal_df.sort_values("gene_symbol").reset_index(drop=True)
            optimal_df.to_csv(optimal_path, index=False)
            self.logger.info(
                "Best group: size=%d, r=%.4f", best_size, best_r,
            )
        else:
            pd.DataFrame(
                columns=["gene_symbol", "log_tpm", "log_counts", "correlation", "group_size"]
            ).to_csv(optimal_path, index=False)
            self.logger.warning(
                "No groups met the correlation threshold (%.4f).",
                cfg.correlation_threshold,
            )
            best_size = 0
            best_r = float("nan")
        output_files.append(str(optimal_path))

        # 4. Detailed Excel results
        xlsx_path = output_dir / "detailed_results.xlsx"
        _save_detailed_results(results_by_size, df, xlsx_path)
        output_files.append(str(xlsx_path))

        # Build result
        metadata = {
            "n_genes_available": n_genes,
            "sizes_tested": len(correlation_trend),
            "sizes_above_threshold": len(results_by_size),
            "best_group_size": best_size,
            "best_correlation": round(best_r, 4) if not np.isnan(best_r) else None,
            "correlation_threshold": cfg.correlation_threshold,
            "merged_counts_file": str(merged_path),
        }

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata=metadata,
        )

        self.write_run_metadata(result, start_time, parameters=metadata)
        return result

    # ------------------------------------------------------------------
    # Path resolution helpers
    # ------------------------------------------------------------------

    def _resolve_merged_counts_file(self) -> Path | None:
        """Resolve the merged counts CSV from the correlation stage output.

        Priority:
        1. If ``distance_threshold`` is set, look for that specific file.
        2. Otherwise, read ``info_{xp}.csv`` and pick the threshold with the
           highest Pearson r.
        3. Fallback: pick the first CSV in ``merged_counts/``.
        """
        output_dir = Path(self.config.paths.output_dir)
        xp_name = self.config.experiment.name
        merged_dir = output_dir / "correlation" / "merged_counts"
        cfg = self.config.optimize_correlation

        if not merged_dir.is_dir():
            return None

        # 1. Explicit threshold
        if cfg.distance_threshold is not None:
            return merged_dir / f"{xp_name}_{cfg.distance_threshold}.csv"

        # 2. Best threshold from info CSV
        info_path = output_dir / "correlation" / f"info_{xp_name}.csv"
        if info_path.exists():
            try:
                info_df = pd.read_csv(info_path)
                if "pearson_r" in info_df.columns and len(info_df) > 0:
                    valid = info_df.dropna(subset=["pearson_r"])
                    if len(valid) > 0:
                        best_idx = valid["pearson_r"].idxmax()
                        best_thresh = valid.loc[best_idx, "threshold"]
                        candidate = merged_dir / f"{xp_name}_{best_thresh}.csv"
                        if candidate.exists():
                            return candidate
            except Exception:
                pass

        # 3. Fallback: first CSV in merged_counts/
        csvs = sorted(merged_dir.glob("*.csv"))
        if csvs:
            return csvs[0]

        return None
