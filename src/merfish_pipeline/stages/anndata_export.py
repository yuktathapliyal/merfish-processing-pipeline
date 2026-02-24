"""``anndata_export`` stage -- export cell-by-gene count matrix to AnnData h5ad.

Converts the ``cell_assignment`` output into a scanpy-compatible AnnData
object for downstream single-cell analysis.

Algorithm
---------
1. Load barcodes CSV with ``Cell_ID`` column (from ``cell_assignment`` or
   explicit config).
2. Map ``barcode_id`` to ``gene_symbol`` via the codebook.
3. Build a cell x gene count matrix using ``pd.crosstab``.
4. Store spatial coordinates (mean ``global_x`` / ``global_y`` per cell) in
   ``adata.obsm['spatial']``.
5. Store cell metadata in ``adata.obs`` (FOV, n_barcodes, n_genes).
6. Store gene metadata in ``adata.var`` (is_blank flag).
7. Optionally filter cells by minimum barcode count and exclude blank genes.

Outputs
-------
- ``{output_dir}/anndata_export/{experiment}.h5ad``
- ``{output_dir}/anndata_export/cell_gene_matrix.csv``
- ``{output_dir}/anndata_export/cell_metadata.csv``
- ``{output_dir}/anndata_export/run_metadata.json``
"""

from __future__ import annotations

import logging
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)

_BLANK_RE = re.compile(r"^[Bb]lank[-_]?\d+$")


def _is_blank(gene_symbol: str) -> bool:
    return bool(_BLANK_RE.match(str(gene_symbol)))


# ---------------------------------------------------------------------------
# Codebook helpers (same pattern as barcode_qc / correlation)
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
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("anndata_export")
class AnnDataExportStage(PipelineStage):
    """Export cell-by-gene count matrix to AnnData h5ad format."""

    description = "Export cell x gene matrix to AnnData (h5ad) and CSV"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []

        barcodes_path = self._resolve_barcodes_path()
        if barcodes_path is None:
            errors.append(
                "Cannot locate barcodes CSV with Cell_ID column. "
                "Set anndata_export.barcodes_file or run cell_assignment first."
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
        xp_name = self.config.experiment.name
        return (
            (out / f"{xp_name}.h5ad").exists()
            or (out / "cell_gene_matrix.csv").exists()
        )

    def run(self, dry_run: bool = False) -> StageResult:
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        cfg = self.config.anndata_export
        xp_name = self.config.experiment.name

        barcodes_path = self._resolve_barcodes_path()
        codebook_path = self._resolve_codebook_path()

        if barcodes_path is None or codebook_path is None:
            return StageResult(
                status="failed",
                error="Could not resolve barcodes or codebook path.",
            )

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would export AnnData from %s with codebook %s",
                barcodes_path,
                codebook_path,
            )
            return StageResult(status="skipped", metadata={"dry_run": True})

        # ----------------------------------------------------------
        # 1. Load data
        # ----------------------------------------------------------
        self.logger.info("Loading barcodes from %s ...", barcodes_path)
        barcodes = read_sheet(barcodes_path)
        self.logger.info("Loaded %d barcode rows.", len(barcodes))

        if "Cell_ID" not in barcodes.columns:
            return StageResult(
                status="failed",
                error=(
                    f"Barcodes file {barcodes_path} does not contain a "
                    "'Cell_ID' column. Run the cell_assignment stage first."
                ),
            )

        codebook = _load_codebook(codebook_path)

        # ----------------------------------------------------------
        # 2. Filter to assigned barcodes (drop background)
        # ----------------------------------------------------------
        assigned = barcodes[barcodes["Cell_ID"].notna()].copy()
        self.logger.info(
            "%d / %d barcodes assigned to cells.",
            len(assigned),
            len(barcodes),
        )

        if assigned.empty:
            return StageResult(
                status="failed",
                error="No barcodes are assigned to cells (all Cell_ID are null).",
            )

        # ----------------------------------------------------------
        # 3. Map barcode_id -> gene_symbol
        # ----------------------------------------------------------
        if "barcode_id" not in assigned.columns:
            return StageResult(
                status="failed",
                error="Barcodes CSV missing 'barcode_id' column.",
            )

        cb_map = codebook.set_index("barcode_id")["gene_symbol"]
        assigned["gene_symbol"] = assigned["barcode_id"].map(cb_map)

        unmapped = assigned["gene_symbol"].isna().sum()
        if unmapped > 0:
            self.logger.warning(
                "%d barcodes could not be mapped to a gene symbol (dropped).",
                unmapped,
            )
            assigned = assigned[assigned["gene_symbol"].notna()]

        # ----------------------------------------------------------
        # 4. Optional: filter by min_barcodes_per_cell
        # ----------------------------------------------------------
        if cfg.min_barcodes_per_cell > 0:
            cell_counts = assigned.groupby("Cell_ID").size()
            keep_cells = cell_counts[
                cell_counts >= cfg.min_barcodes_per_cell
            ].index
            n_before = assigned["Cell_ID"].nunique()
            assigned = assigned[assigned["Cell_ID"].isin(keep_cells)]
            n_after = assigned["Cell_ID"].nunique()
            self.logger.info(
                "Filtered cells with < %d barcodes: %d -> %d cells.",
                cfg.min_barcodes_per_cell,
                n_before,
                n_after,
            )

        # ----------------------------------------------------------
        # 5. Build count matrix
        # ----------------------------------------------------------
        count_matrix = pd.crosstab(
            assigned["Cell_ID"], assigned["gene_symbol"]
        )
        self.logger.info(
            "Count matrix: %d cells x %d genes.",
            count_matrix.shape[0],
            count_matrix.shape[1],
        )

        # ----------------------------------------------------------
        # 6. Build cell metadata (obs)
        # ----------------------------------------------------------
        cell_meta = assigned.groupby("Cell_ID").agg(
            n_barcodes=("gene_symbol", "size"),
            n_genes=("gene_symbol", "nunique"),
        )

        # Detect FOV column
        fov_col = None
        for c in ["fov", "FOV", "Fov"]:
            if c in assigned.columns:
                fov_col = c
                break

        if fov_col is not None:
            # Most common FOV per cell
            fov_mode = (
                assigned.groupby("Cell_ID")[fov_col]
                .agg(lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else x.iloc[0])
            )
            cell_meta["fov"] = fov_mode

        # ----------------------------------------------------------
        # 7. Spatial coordinates (mean global_x, global_y per cell)
        # ----------------------------------------------------------
        spatial_coords = None
        if "global_x" in assigned.columns and "global_y" in assigned.columns:
            spatial_coords = (
                assigned.groupby("Cell_ID")[["global_x", "global_y"]]
                .mean()
            )

        # ----------------------------------------------------------
        # 8. Gene metadata (var)
        # ----------------------------------------------------------
        gene_meta = pd.DataFrame(
            {"is_blank": [_is_blank(g) for g in count_matrix.columns]},
            index=count_matrix.columns,
        )

        # ----------------------------------------------------------
        # 9. Optionally exclude blanks from count matrix
        # ----------------------------------------------------------
        if cfg.exclude_blanks:
            blank_genes = gene_meta.index[gene_meta["is_blank"]]
            n_blank = len(blank_genes)
            if n_blank > 0:
                count_matrix = count_matrix.drop(columns=blank_genes)
                gene_meta = gene_meta.loc[~gene_meta["is_blank"]]
                self.logger.info("Excluded %d blank genes from matrix.", n_blank)

        # ----------------------------------------------------------
        # 10. Write outputs
        # ----------------------------------------------------------
        output_dir.mkdir(parents=True, exist_ok=True)
        output_files: list[str] = []

        # CSV: cell x gene count matrix
        matrix_path = output_dir / "cell_gene_matrix.csv"
        count_matrix.to_csv(matrix_path)
        output_files.append(str(matrix_path))
        self.logger.info("Wrote count matrix: %s", matrix_path)

        # CSV: cell metadata
        meta_path = output_dir / "cell_metadata.csv"
        cell_meta.to_csv(meta_path)
        output_files.append(str(meta_path))

        # h5ad: AnnData (optional -- requires anndata package)
        h5ad_path = output_dir / f"{xp_name}.h5ad"
        try:
            import anndata

            adata = anndata.AnnData(
                X=count_matrix.values,
                obs=cell_meta.reindex(count_matrix.index),
                var=gene_meta,
            )
            adata.obs_names = count_matrix.index.astype(str)
            adata.var_names = count_matrix.columns.astype(str)

            if spatial_coords is not None:
                spatial_arr = (
                    spatial_coords.reindex(count_matrix.index)
                    .values.astype(np.float32)
                )
                adata.obsm["spatial"] = spatial_arr

            adata.write_h5ad(h5ad_path)
            output_files.append(str(h5ad_path))
            self.logger.info("Wrote AnnData: %s", h5ad_path)
        except ImportError:
            self.logger.warning(
                "anndata package not installed -- skipping h5ad export. "
                "Install with: pip install 'merfish-pipeline[export]'"
            )
        except Exception as exc:
            self.logger.error("Failed to write h5ad: %s", exc)

        # ----------------------------------------------------------
        # 11. Summary
        # ----------------------------------------------------------
        metadata = {
            "n_cells": count_matrix.shape[0],
            "n_genes": count_matrix.shape[1],
            "n_barcodes_used": len(assigned),
            "exclude_blanks": cfg.exclude_blanks,
            "min_barcodes_per_cell": cfg.min_barcodes_per_cell,
            "has_spatial": spatial_coords is not None,
            "h5ad_written": h5ad_path.exists(),
        }

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata=metadata,
        )

        self.write_run_metadata(
            result,
            start_time,
            parameters={
                "barcodes_file": str(barcodes_path),
                "codebook_file": str(codebook_path),
                **metadata,
            },
        )

        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_barcodes_path(self) -> Path | None:
        """Resolve barcodes CSV: explicit -> cell_assignment -> filter_barcodes."""
        explicit = self.config.anndata_export.barcodes_file
        if explicit is not None:
            return Path(explicit)

        output_dir = Path(self.config.paths.output_dir)

        # cell_assignment output (preferred — has Cell_ID)
        ca_output = output_dir / "cell_assignment" / "barcodes_assigned.csv"
        if ca_output.exists():
            return ca_output

        # filter_barcodes (may not have Cell_ID — stage will fail gracefully)
        fb_output = output_dir / "filter_barcodes" / "barcodes_filtered.csv"
        if fb_output.exists():
            return fb_output

        return None

    def _resolve_codebook_path(self) -> Path | None:
        if self.config.merlin.codebook_template is not None:
            return Path(self.config.merlin.codebook_template)
        return None
