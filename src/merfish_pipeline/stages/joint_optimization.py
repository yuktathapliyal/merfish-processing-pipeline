"""``joint_optimization`` stage -- find gene subgroups maximizing correlation
across multiple experiments simultaneously.

Uses simulated annealing to find, for each target group size, the subset of
genes that maximizes the *average* Pearson correlation between merFISH barcode
log-counts and bulk RNA-seq log-TPM across all experiments.

The current experiment is always included.  Additional experiments are
referenced via the ``joint_experiments`` config list, which points to their
``correlation`` stage output directories (or directly to merged-counts CSVs).

Outputs
-------
- ``{output_dir}/joint_optimization/correlation_trend.csv``
- ``{output_dir}/joint_optimization/correlation_trend.png``
- ``{output_dir}/joint_optimization/optimal_genes.csv``
- ``{output_dir}/joint_optimization/detailed_results.xlsx``
- ``{output_dir}/joint_optimization/subgroup_correlations_unlabeled.pdf``
- ``{output_dir}/joint_optimization/subgroup_correlations_labeled.pdf``
- ``{output_dir}/joint_optimization/run_metadata.json``
"""

from __future__ import annotations

import logging
import math
import random
import time
from datetime import datetime
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.optimize_correlation import (
    _calculate_correlation,
    _plot_subgroup_scatter,
)
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Merged-counts resolution helpers
# ---------------------------------------------------------------------------


def _resolve_own_merged_counts(
    output_dir: Path,
    xp_name: str,
    distance_threshold: float | None,
) -> Path | None:
    """Resolve the merged counts CSV for the current experiment.

    Same algorithm as ``OptimizeCorrelationStage._resolve_merged_counts_file``
    but as a standalone function so it can be called from any stage.
    """
    merged_dir = output_dir / "correlation" / "merged_counts"
    if not merged_dir.is_dir():
        return None

    # 1. Explicit threshold
    if distance_threshold is not None:
        return merged_dir / f"{xp_name}_{distance_threshold}.csv"

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
    return csvs[0] if csvs else None


def _resolve_external_merged_counts(
    entry: object,
) -> tuple[str, Path]:
    """Resolve experiment name and merged-counts CSV for an external experiment.

    Parameters
    ----------
    entry
        A ``JointExperimentEntry`` config object with ``name``,
        ``correlation_dir``, ``merged_counts``, and ``distance_threshold``.

    Returns
    -------
    (experiment_name, csv_path)

    Raises
    ------
    FileNotFoundError
        If the CSV cannot be located — message is user-actionable.
    """
    # Priority 1: direct merged_counts path
    if entry.merged_counts is not None:
        p = Path(entry.merged_counts)
        if not p.exists():
            raise FileNotFoundError(
                f"Merged counts file not found: {p}. "
                f"Check that the file exists."
            )
        name = entry.name or p.stem.rsplit("_", 1)[0]
        return name, p

    # Priority 2: auto-discover from correlation_dir
    corr_dir = Path(entry.correlation_dir)  # type: ignore[arg-type]
    if not corr_dir.is_dir():
        raise FileNotFoundError(
            f"Correlation directory not found: {corr_dir}. "
            f"Check that the path exists and the correlation stage has been "
            f"run for that experiment."
        )

    # Find the info CSV
    info_csvs = sorted(corr_dir.glob("info_*.csv"))
    if len(info_csvs) == 0:
        raise FileNotFoundError(
            f"No info_*.csv found in {corr_dir}. "
            f"Run the correlation stage for that experiment first, or "
            f"provide 'merged_counts' directly."
        )
    if len(info_csvs) > 1:
        names = [p.name for p in info_csvs]
        raise FileNotFoundError(
            f"Multiple info CSVs found in {corr_dir}: {names}. "
            f"Provide 'merged_counts' directly to select the correct one."
        )

    info_path = info_csvs[0]
    xp_name = info_path.stem.removeprefix("info_")
    merged_dir = corr_dir / "merged_counts"

    # Explicit threshold override
    if entry.distance_threshold is not None:
        candidate = merged_dir / f"{xp_name}_{entry.distance_threshold}.csv"
        if not candidate.exists():
            raise FileNotFoundError(
                f"Merged counts file not found: {candidate}. "
                f"Available files: {sorted(p.name for p in merged_dir.glob('*.csv'))}"
            )
        name = entry.name or xp_name
        return name, candidate

    # Auto-select best threshold from info CSV
    try:
        info_df = pd.read_csv(info_path)
        if "pearson_r" in info_df.columns and len(info_df) > 0:
            valid = info_df.dropna(subset=["pearson_r"])
            if len(valid) > 0:
                best_idx = valid["pearson_r"].idxmax()
                best_thresh = valid.loc[best_idx, "threshold"]
                candidate = merged_dir / f"{xp_name}_{best_thresh}.csv"
                if candidate.exists():
                    name = entry.name or xp_name
                    return name, candidate
            else:
                raise FileNotFoundError(
                    f"Experiment '{xp_name}' info CSV has no valid Pearson "
                    f"correlations at any threshold. Set "
                    f"'distance_threshold' explicitly in the joint_experiments "
                    f"entry for this experiment."
                )
    except FileNotFoundError:
        raise
    except Exception:
        pass

    # Fallback: first CSV in merged_counts/
    if merged_dir.is_dir():
        csvs = sorted(merged_dir.glob("*.csv"))
        if csvs:
            name = entry.name or xp_name
            return name, csvs[0]

    raise FileNotFoundError(
        f"Could not locate merged counts in {corr_dir}. "
        f"Provide 'merged_counts' directly."
    )


# ---------------------------------------------------------------------------
# Gene alignment
# ---------------------------------------------------------------------------


def _align_experiments(
    dfs: dict[str, pd.DataFrame],
) -> tuple[dict[str, pd.DataFrame], list[str]]:
    """Find common genes across all experiments and align dataframes.

    Returns ``(aligned_dfs, common_gene_list)`` where every dataframe has
    rows in the same order (sorted alphabetically by ``gene_symbol``).
    """
    gene_sets = {
        name: set(df["gene_symbol"].unique()) for name, df in dfs.items()
    }
    common = set.intersection(*gene_sets.values())
    common_genes = sorted(common)

    for name, gs in gene_sets.items():
        n_total = len(gs)
        n_dropped = n_total - len(common)
        logger.info(
            "  %s: %d genes total, %d in common set (%d dropped)",
            name, n_total, len(common), n_dropped,
        )

    if len(common_genes) == 0:
        exp_summary = ", ".join(
            f"{name} ({len(gs)} genes)" for name, gs in gene_sets.items()
        )
        raise ValueError(
            f"No common genes found across experiments: {exp_summary}. "
            f"This usually means the experiments used different codebooks."
        )

    aligned: dict[str, pd.DataFrame] = {}
    for name, df in dfs.items():
        sub = df[df["gene_symbol"].isin(common)].copy()
        sub = sub.drop_duplicates(subset=["gene_symbol"])
        sub = sub.set_index("gene_symbol").loc[common_genes].reset_index()
        aligned[name] = sub

    return aligned, common_genes


# ---------------------------------------------------------------------------
# Multi-experiment simulated annealing
# ---------------------------------------------------------------------------


def _calculate_multi_correlation(
    experiments: list[tuple[np.ndarray, np.ndarray]],
    indices: set[int],
) -> float:
    """Average Pearson *r* across all experiments for the gene subset.

    Returns ``-inf`` if **any** experiment fails (ensures the subset works
    for all experiments).
    """
    rs: list[float] = []
    for log_tpm, log_counts in experiments:
        r = _calculate_correlation(log_tpm, log_counts, indices)
        if r == float("-inf"):
            return float("-inf")
        rs.append(r)
    return sum(rs) / len(rs)


def _optimize_joint_group(
    experiments: list[tuple[np.ndarray, np.ndarray]],
    n_genes: int,
    target_size: int,
    *,
    max_iterations: int,
    initial_temperature: float,
    cooling_rate: float,
    rng: random.Random,
) -> tuple[set[int], float]:
    """Single SA run for a given target group size using multi-experiment scoring.

    Returns ``(best_indices, best_avg_correlation)``.
    """
    all_indices = list(range(n_genes))
    current = set(rng.sample(all_indices, min(target_size, n_genes)))
    current_r = _calculate_multi_correlation(experiments, current)

    best = set(current)
    best_r = current_r
    temperature = initial_temperature

    for _ in range(max_iterations):
        candidate = set(current)
        outside = [i for i in all_indices if i not in candidate]
        if not outside:
            break

        remove_gene = rng.choice(list(candidate))
        add_gene = rng.choice(outside)
        candidate.discard(remove_gene)
        candidate.add(add_gene)

        candidate_r = _calculate_multi_correlation(experiments, candidate)
        delta = candidate_r - current_r

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


def _find_progressive_joint_groups(
    experiments: list[tuple[np.ndarray, np.ndarray]],
    experiment_names: list[str],
    n_genes: int,
    *,
    size_range: range,
    correlation_threshold: float,
    n_attempts: int,
    max_iterations: int,
    initial_temperature: float,
    cooling_rate: float,
    rng: random.Random,
) -> tuple[
    dict[int, tuple[set[int], float, dict[str, float]]],
    list[tuple[int, float, dict[str, float]]],
    dict[int, tuple[set[int], float, dict[str, float]]],
]:
    """Run joint SA optimization across a range of target sizes.

    Returns
    -------
    results_by_size : dict
        ``{size: (best_indices, best_avg_r, {xp: r})}`` for sizes where
        the threshold was met.
    correlation_trend : list
        ``[(size, best_avg_r, {xp: r}), ...]`` for all tested sizes.
    all_results_by_size : dict
        ``{size: (best_indices, best_avg_r, {xp: r})}`` for every tested
        size regardless of threshold.
    """
    results_by_size: dict[int, tuple[set[int], float, dict[str, float]]] = {}
    correlation_trend: list[tuple[int, float, dict[str, float]]] = []
    all_results_by_size: dict[int, tuple[set[int], float, dict[str, float]]] = {}

    total_sizes = len(list(size_range))
    t0 = time.monotonic()

    for size_idx, target_size in enumerate(size_range, 1):
        if target_size > n_genes:
            logger.info(
                "  Skipping size %d (exceeds %d common genes)",
                target_size, n_genes,
            )
            break

        best_indices: set[int] | None = None
        best_avg_r = float("-inf")

        for _attempt in range(n_attempts):
            indices, avg_r = _optimize_joint_group(
                experiments, n_genes, target_size,
                max_iterations=max_iterations,
                initial_temperature=initial_temperature,
                cooling_rate=cooling_rate,
                rng=rng,
            )
            if avg_r > best_avg_r:
                best_avg_r = avg_r
                best_indices = indices

        # Compute per-experiment correlations for the best subset
        per_exp: dict[str, float] = {}
        if best_indices is not None:
            idx_list = list(best_indices)
            for name, (log_tpm, log_counts) in zip(
                experiment_names, experiments
            ):
                try:
                    r, _ = pearsonr(log_tpm[idx_list], log_counts[idx_list])
                    per_exp[name] = float(r) if not np.isnan(r) else float("nan")
                except Exception:
                    per_exp[name] = float("nan")

        elapsed = time.monotonic() - t0
        logger.info(
            "  Size %3d (%d/%d): best avg r = %.4f across %d experiments "
            "[elapsed: %.0fs]",
            target_size, size_idx, total_sizes, best_avg_r,
            len(experiments), elapsed,
        )

        correlation_trend.append((target_size, best_avg_r, per_exp))

        if best_indices is not None:
            all_results_by_size[target_size] = (best_indices, best_avg_r, per_exp)

        if best_avg_r >= correlation_threshold and best_indices is not None:
            results_by_size[target_size] = (best_indices, best_avg_r, per_exp)

    return results_by_size, correlation_trend, all_results_by_size


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def _plot_joint_correlation_trend(
    trend_df: pd.DataFrame,
    experiment_names: list[str],
    threshold: float,
    output_path: Path,
) -> None:
    """Line plot with one line per experiment, a bold average line, and threshold."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Per-experiment lines
    colors = plt.cm.tab10.colors  # type: ignore[attr-defined]
    for i, name in enumerate(experiment_names):
        col = f"{name}_r"
        if col in trend_df.columns:
            color = colors[i % len(colors)]
            ax.plot(
                trend_df["group_size"], trend_df[col],
                "o-", color=color, linewidth=1, markersize=4,
                alpha=0.7, label=name,
            )

    # Average line (bold)
    ax.plot(
        trend_df["group_size"], trend_df["avg_correlation"],
        "s-", color="black", linewidth=2.5, markersize=5,
        label="Average",
    )

    # Threshold
    ax.axhline(
        threshold, color="red", linestyle="--", alpha=0.7,
        label=f"Threshold = {threshold}",
    )

    ax.set_xlabel("Gene Group Size", fontsize=12)
    ax.set_ylabel("Pearson Correlation", fontsize=12)
    ax.set_title(
        "Joint Correlation vs Gene Group Size (Simulated Annealing)",
        fontsize=13,
    )
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_joint_detailed_results(
    all_results_by_size: dict[int, tuple[set[int], float, dict[str, float]]],
    aligned_dfs: dict[str, pd.DataFrame],
    common_genes: list[str],
    experiment_names: list[str],
    output_path: Path,
) -> None:
    """Write Excel workbook with summary sheet and per-size gene sheets."""
    summary_rows: list[dict] = []
    sheets: dict[str, pd.DataFrame] = {}

    for size in sorted(all_results_by_size):
        indices, avg_r, per_exp = all_results_by_size[size]
        idx_list = sorted(indices)
        gene_names = [common_genes[i] for i in idx_list]

        row: dict = {
            "group_size": size,
            "avg_correlation": round(avg_r, 4),
            "n_genes": len(indices),
        }
        for name in experiment_names:
            row[f"{name}_r"] = round(per_exp.get(name, float("nan")), 4)

        # Total counts per experiment for this subset
        for name in experiment_names:
            df = aligned_dfs[name]
            subset = df[df["gene_symbol"].isin(gene_names)]
            if "counts" in subset.columns:
                row[f"{name}_total_counts"] = int(subset["counts"].sum())

        row["gene_list"] = ", ".join(sorted(gene_names))
        summary_rows.append(row)

        # Per-size gene sheet
        gene_data = pd.DataFrame({"gene_symbol": sorted(gene_names)})
        for name in experiment_names:
            df = aligned_dfs[name]
            subset = df[df["gene_symbol"].isin(gene_names)].copy()
            subset = subset.set_index("gene_symbol")
            gene_data[f"{name}_log_tpm"] = (
                gene_data["gene_symbol"].map(subset["log_tpm"])
            )
            gene_data[f"{name}_log_counts"] = (
                gene_data["gene_symbol"].map(subset["log_counts"])
            )

        sheet_name = f"Group_{size}"
        if len(sheet_name) > 31:
            sheet_name = sheet_name[:31]
        sheets[sheet_name] = gene_data

    if not summary_rows:
        summary_rows.append({
            "group_size": 0,
            "avg_correlation": float("nan"),
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


@register_stage("joint_optimization")
class JointOptimizationStage(PipelineStage):
    """Find gene subgroups maximizing bulk correlation across multiple experiments."""

    description = (
        "Joint simulated annealing optimization of gene subgroups "
        "for maximum correlation with bulk RNA-seq across experiments"
    )

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []
        cfg = self.config.joint_optimization

        if not cfg.joint_experiments:
            errors.append(
                "joint_optimization is enabled but no experiments are listed "
                "in 'joint_experiments'. Add at least one experiment entry."
            )

        # Validate each external experiment entry exists on disk
        for i, entry in enumerate(cfg.joint_experiments):
            if entry.merged_counts is not None:
                if not Path(entry.merged_counts).exists():
                    errors.append(
                        f"joint_experiments[{i}]: merged counts file not "
                        f"found: {entry.merged_counts}"
                    )
            elif entry.correlation_dir is not None:
                if not Path(entry.correlation_dir).is_dir():
                    errors.append(
                        f"joint_experiments[{i}]: correlation directory not "
                        f"found: {entry.correlation_dir}. Check that the "
                        f"path exists and the correlation stage has been run "
                        f"for that experiment."
                    )

        # Validate current experiment's merged counts
        output_dir = Path(self.config.paths.output_dir)
        xp_name = self.config.experiment.name
        own_path = _resolve_own_merged_counts(
            output_dir, xp_name, cfg.distance_threshold,
        )
        if own_path is None:
            errors.append(
                "Cannot locate merged counts for the current experiment. "
                "Run the 'correlation' stage first, or set "
                "joint_optimization.distance_threshold explicitly."
            )
        elif not own_path.exists():
            errors.append(
                f"Merged counts file not found for current experiment: "
                f"{own_path}"
            )

        # SA parameter validation
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
            and (out / "subgroup_correlations_unlabeled.pdf").exists()
            and (out / "subgroup_correlations_labeled.pdf").exists()
        )

    def run(self, dry_run: bool = False) -> StageResult:
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        cfg = self.config.joint_optimization
        base_output_dir = Path(self.config.paths.output_dir)
        xp_name = self.config.experiment.name

        # --- Resolve all merged-counts CSVs ---

        own_path = _resolve_own_merged_counts(
            base_output_dir, xp_name, cfg.distance_threshold,
        )
        if own_path is None or not own_path.exists():
            return StageResult(
                status="failed",
                error="Could not resolve merged counts for current experiment.",
            )

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would run joint optimization with %d external "
                "experiments",
                len(cfg.joint_experiments),
            )
            return StageResult(status="skipped", metadata={"dry_run": True})

        # Resolve external experiments
        experiment_csvs: dict[str, Path] = {xp_name: own_path}
        for i, entry in enumerate(cfg.joint_experiments):
            try:
                name, csv_path = _resolve_external_merged_counts(entry)
            except FileNotFoundError as exc:
                return StageResult(
                    status="failed",
                    error=f"joint_experiments[{i}]: {exc}",
                )
            if name in experiment_csvs:
                return StageResult(
                    status="failed",
                    error=(
                        f"Duplicate experiment name '{name}' resolved from "
                        f"two different paths. Add explicit 'name' fields "
                        f"to the joint_experiments entries to distinguish them."
                    ),
                )
            experiment_csvs[name] = csv_path

        # --- Load all CSVs ---

        self.logger.info(
            "Loading merged counts for %d experiments...",
            len(experiment_csvs),
        )
        raw_dfs: dict[str, pd.DataFrame] = {}
        for name, csv_path in experiment_csvs.items():
            self.logger.info("  %s: %s", name, csv_path)
            df = read_sheet(csv_path)
            for col in ("log_tpm", "log_counts", "gene_symbol"):
                if col not in df.columns:
                    return StageResult(
                        status="failed",
                        error=(
                            f"Experiment '{name}' merged counts missing "
                            f"column '{col}'. Available: {list(df.columns)}"
                        ),
                    )
            df = df.dropna(
                subset=["log_tpm", "log_counts", "gene_symbol"]
            ).reset_index(drop=True)
            raw_dfs[name] = df

        # --- Align to common gene set ---

        self.logger.info("Aligning experiments to common gene set...")
        try:
            aligned_dfs, common_genes = _align_experiments(raw_dfs)
        except ValueError as exc:
            return StageResult(status="failed", error=str(exc))

        n_common = len(common_genes)
        self.logger.info("Common gene set: %d genes", n_common)

        # Cap size range at available genes
        effective_end = min(cfg.size_range_end, n_common) + 1
        size_range = range(cfg.size_range_start, effective_end, cfg.size_range_step)
        if not list(size_range):
            return StageResult(
                status="failed",
                error=(
                    f"Empty size range: start={cfg.size_range_start}, "
                    f"end={effective_end - 1}, step={cfg.size_range_step}, "
                    f"common genes={n_common}"
                ),
            )

        if n_common < cfg.size_range_start:
            self.logger.warning(
                "Only %d common genes — fewer than size_range_start (%d). "
                "Joint optimization may not be meaningful.",
                n_common, cfg.size_range_start,
            )

        # --- Build experiment arrays ---

        experiment_names = list(aligned_dfs.keys())
        experiments: list[tuple[np.ndarray, np.ndarray]] = []
        for name in experiment_names:
            df = aligned_dfs[name]
            experiments.append((df["log_tpm"].values, df["log_counts"].values))

        # --- Run joint SA ---

        rng = random.Random(cfg.random_seed)
        self.logger.info(
            "Running joint SA optimization: sizes %d-%d (step %d), "
            "%d attempts, %d iterations each, %d experiments",
            cfg.size_range_start, effective_end - 1, cfg.size_range_step,
            cfg.n_attempts, cfg.max_iterations, len(experiments),
        )

        results_by_size, correlation_trend, all_results_by_size = (
            _find_progressive_joint_groups(
                experiments, experiment_names, n_common,
                size_range=size_range,
                correlation_threshold=cfg.correlation_threshold,
                n_attempts=cfg.n_attempts,
                max_iterations=cfg.max_iterations,
                initial_temperature=cfg.initial_temperature,
                cooling_rate=cfg.cooling_rate,
                rng=rng,
            )
        )

        # --- Write outputs ---

        output_dir.mkdir(parents=True, exist_ok=True)
        output_files: list[str] = []

        # 1. Correlation trend CSV
        trend_rows: list[dict] = []
        for size, avg_r, per_exp in correlation_trend:
            row: dict = {
                "group_size": size,
                "avg_correlation": round(avg_r, 4),
            }
            for name in experiment_names:
                row[f"{name}_r"] = round(per_exp.get(name, float("nan")), 4)
            # Total counts per experiment at this size
            if size in all_results_by_size:
                indices = all_results_by_size[size][0]
                gene_names = [common_genes[i] for i in sorted(indices)]
                for name in experiment_names:
                    df = aligned_dfs[name]
                    subset = df[df["gene_symbol"].isin(gene_names)]
                    if "counts" in subset.columns:
                        row[f"{name}_total_counts"] = int(subset["counts"].sum())
            trend_rows.append(row)

        trend_df = pd.DataFrame(trend_rows)
        trend_path = output_dir / "correlation_trend.csv"
        trend_df.to_csv(trend_path, index=False)
        output_files.append(str(trend_path))

        # 2. Correlation trend plot
        trend_plot_path = output_dir / "correlation_trend.png"
        _plot_joint_correlation_trend(
            trend_df, experiment_names,
            cfg.correlation_threshold, trend_plot_path,
        )
        output_files.append(str(trend_plot_path))

        # 3. Optimal genes CSV
        optimal_path = output_dir / "optimal_genes.csv"
        if results_by_size:
            best_size = max(
                results_by_size,
                key=lambda s: results_by_size[s][1],
            )
            best_indices, best_avg_r, best_per_exp = results_by_size[best_size]
            gene_names = [common_genes[i] for i in sorted(best_indices)]
            opt_rows: list[dict] = []
            for gene in gene_names:
                opt_row: dict = {"gene_symbol": gene}
                for name in experiment_names:
                    df = aligned_dfs[name]
                    match = df[df["gene_symbol"] == gene]
                    if len(match) > 0:
                        opt_row[f"{name}_log_tpm"] = match.iloc[0]["log_tpm"]
                        opt_row[f"{name}_log_counts"] = match.iloc[0]["log_counts"]
                opt_rows.append(opt_row)
            optimal_df = pd.DataFrame(opt_rows)
            optimal_df["avg_correlation"] = best_avg_r
            for name in experiment_names:
                optimal_df[f"{name}_r"] = best_per_exp.get(name, float("nan"))
            optimal_df["group_size"] = best_size
            optimal_df.to_csv(optimal_path, index=False)
            self.logger.info(
                "Best joint group: size=%d, avg r=%.4f", best_size, best_avg_r,
            )
        else:
            pd.DataFrame(columns=["gene_symbol", "avg_correlation", "group_size"]).to_csv(
                optimal_path, index=False,
            )
            self.logger.warning(
                "No groups met the correlation threshold (%.4f).",
                cfg.correlation_threshold,
            )
            best_size = 0
            best_avg_r = float("nan")
        output_files.append(str(optimal_path))

        # 4. Detailed Excel results
        xlsx_path = output_dir / "detailed_results.xlsx"
        _save_joint_detailed_results(
            all_results_by_size, aligned_dfs, common_genes,
            experiment_names, xlsx_path,
        )
        output_files.append(str(xlsx_path))

        # 5. Per-subgroup scatter PDFs (multi-panel: one subplot per experiment)
        n_exp = len(experiment_names)
        unlabeled_pdf_path = output_dir / "subgroup_correlations_unlabeled.pdf"
        labeled_pdf_path = output_dir / "subgroup_correlations_labeled.pdf"

        with PdfPages(str(unlabeled_pdf_path)) as pdf_u, \
             PdfPages(str(labeled_pdf_path)) as pdf_l:
            for size in sorted(all_results_by_size):
                indices, avg_r, per_exp = all_results_by_size[size]
                gene_names = [common_genes[i] for i in sorted(indices)]

                for pdf_obj, with_labels in ((pdf_u, False), (pdf_l, True)):
                    fig, axes = plt.subplots(
                        1, n_exp,
                        figsize=(9 * n_exp, 9),
                        squeeze=False,
                    )
                    for j, name in enumerate(experiment_names):
                        ax = axes[0, j]
                        df = aligned_dfs[name]
                        sub_df = df[df["gene_symbol"].isin(gene_names)][
                            ["gene_symbol", "log_tpm", "log_counts"]
                        ].copy()
                        exp_r = per_exp.get(name, float("nan"))
                        _plot_subgroup_scatter(
                            sub_df, ax,
                            group_size=size,
                            best_r=exp_r,
                            xp_name=name,
                            with_labels=with_labels,
                        )
                    fig.suptitle(
                        f"Joint subgroup size {size} "
                        f"(avg r = {avg_r:.3f})",
                        fontsize=14, y=1.01,
                    )
                    fig.tight_layout()
                    pdf_obj.savefig(fig)
                    plt.close(fig)

        output_files.append(str(unlabeled_pdf_path))
        output_files.append(str(labeled_pdf_path))
        self.logger.info(
            "Wrote per-subgroup scatter PDFs (%d pages each, %d experiments)",
            len(all_results_by_size), n_exp,
        )

        # --- Build result ---

        metadata = {
            "experiments": experiment_names,
            "n_experiments": len(experiment_names),
            "n_common_genes": n_common,
            "sizes_tested": len(correlation_trend),
            "sizes_above_threshold": len(results_by_size),
            "best_group_size": best_size,
            "best_avg_correlation": (
                round(best_avg_r, 4) if not np.isnan(best_avg_r) else None
            ),
            "correlation_threshold": cfg.correlation_threshold,
        }

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata=metadata,
        )
        self.write_run_metadata(result, start_time, parameters=metadata)
        return result
