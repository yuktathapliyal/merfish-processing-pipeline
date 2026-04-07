"""``spatial_visualization`` stage -- interactive 3D barcode scatter plot.

Produces a single interactive plotly HTML file with per-FOV 3D scatter plots
of decoded barcodes.  Two views: gene-colored (barcodes colored by gene
identity) and cell-colored (barcodes colored by Cell_ID, unassigned dimmed).
A dropdown selector switches between FOVs and views.

Outputs
-------
- ``{output_dir}/spatial_visualization/spatial_3d.html``
- ``{output_dir}/spatial_visualization/run_metadata.json``
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

logger = logging.getLogger(__name__)

# Column auto-detection candidates
_GLOBAL_X_CANDIDATES = ["global_x", "x", "X"]
_GLOBAL_Y_CANDIDATES = ["global_y", "y", "Y"]
_Z_CANDIDATES = ["z", "Z", "zIndex", "z_index", "zPos", "zpos"]
_FOV_CANDIDATES = ["fov", "FOV", "Fov"]


def _detect_column(df: pd.DataFrame, candidates: list[str], label: str) -> str:
    """Return the first matching column name from *candidates*."""
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(
        f"Cannot auto-detect {label} column. "
        f"Tried {candidates}; available: {list(df.columns)}"
    )


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


def _build_figure(
    barcodes: pd.DataFrame,
    x_col: str,
    y_col: str,
    z_col: str,
    fov_col: str,
    has_cell_id: bool,
    cfg,
):
    """Build the plotly figure with gene and cell views per FOV."""
    import plotly.graph_objects as go
    import plotly.colors

    fovs = sorted(barcodes[fov_col].unique())

    # Gene color mapping
    unique_genes = sorted(barcodes["gene_symbol"].dropna().unique())
    palette = plotly.colors.qualitative.Alphabet  # 26 colors
    gene_color_map = {g: palette[i % len(palette)] for i, g in enumerate(unique_genes)}

    gene_traces = []
    cell_traces = []

    for fov in fovs:
        fov_df = barcodes[barcodes[fov_col] == fov]

        # Gene-colored trace
        colors = fov_df["gene_symbol"].map(gene_color_map).fillna("#888888").tolist()
        hover = [
            f"Gene: {g}<br>FOV: {fov}<br>x: {x:.1f}, y: {y:.1f}, z: {z}"
            for g, x, y, z in zip(
                fov_df["gene_symbol"],
                fov_df[x_col],
                fov_df[y_col],
                fov_df[z_col],
            )
        ]

        gene_traces.append(go.Scatter3d(
            x=fov_df[x_col].values,
            y=fov_df[y_col].values,
            z=fov_df[z_col].values,
            mode="markers",
            marker=dict(size=cfg.marker_size, color=colors, opacity=0.7),
            hovertext=hover,
            hoverinfo="text",
            name=f"Gene FOV {fov}",
            showlegend=False,
        ))

        # Cell-colored trace (if Cell_ID present)
        if has_cell_id:
            cell_colors = []
            cell_sizes = []
            cell_hover = []

            for _, row in fov_df.iterrows():
                cid = row.get("Cell_ID")
                gene = row.get("gene_symbol", "?")
                if pd.notna(cid):
                    cidx = hash(str(cid)) % len(palette)
                    cell_colors.append(palette[cidx])
                    cell_sizes.append(cfg.marker_size)
                    cell_hover.append(
                        f"Cell: {cid}<br>Gene: {gene}<br>FOV: {fov}"
                    )
                else:
                    cell_colors.append(cfg.unassigned_color)
                    cell_sizes.append(cfg.unassigned_marker_size)
                    cell_hover.append(
                        f"Unassigned<br>Gene: {gene}<br>FOV: {fov}"
                    )

            cell_traces.append(go.Scatter3d(
                x=fov_df[x_col].values,
                y=fov_df[y_col].values,
                z=fov_df[z_col].values,
                mode="markers",
                marker=dict(size=cell_sizes, color=cell_colors, opacity=0.7),
                hovertext=cell_hover,
                hoverinfo="text",
                name=f"Cell FOV {fov}",
                showlegend=False,
            ))

    # Combine: [gene_fov0, gene_fov1, ..., cell_fov0, cell_fov1, ...]
    all_traces = gene_traces + cell_traces
    fig = go.Figure(data=all_traces)

    n_gene = len(gene_traces)
    n_cell = len(cell_traces)
    n_total = n_gene + n_cell

    # Build dropdown buttons
    buttons = []

    # Gene view entries
    for i, fov in enumerate(fovs):
        vis = [False] * n_total
        vis[i] = True
        buttons.append(dict(
            label=f"Gene (FOV {fov})",
            method="update",
            args=[{"visible": vis}],
        ))

    # Cell view entries (if available)
    if has_cell_id:
        for i, fov in enumerate(fovs):
            vis = [False] * n_total
            vis[n_gene + i] = True
            buttons.append(dict(
                label=f"Cell (FOV {fov})",
                method="update",
                args=[{"visible": vis}],
            ))

    # Set initial visibility: first FOV, gene view
    initial_vis = [False] * n_total
    if n_gene > 0:
        initial_vis[0] = True
    for trace, v in zip(fig.data, initial_vis):
        trace.visible = v

    fig.update_layout(
        updatemenus=[dict(
            type="dropdown",
            direction="down",
            buttons=buttons,
            x=0.0,
            y=1.15,
            xanchor="left",
            yanchor="top",
            showactive=True,
        )],
        scene=dict(
            xaxis_title=x_col,
            yaxis_title=y_col,
            zaxis_title=z_col,
            bgcolor="black",
        ),
        paper_bgcolor="black",
        font=dict(color="white"),
        title="Spatial Barcode Visualization",
        height=800,
        width=1000,
    )

    return fig


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("spatial_visualization")
class SpatialVisualizationStage(PipelineStage):
    """Interactive 3D spatial barcode visualization (plotly HTML)."""

    description = "Generate interactive 3D scatter plot of decoded barcodes"

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []

        barcodes_path = self._resolve_barcodes_path()
        if barcodes_path is None:
            errors.append(
                "Cannot locate barcodes CSV. Set spatial_visualization.barcodes_file "
                "or ensure cell_assignment / filter_barcodes / MERlin output exists."
            )
        elif not barcodes_path.exists():
            errors.append(f"Barcodes file does not exist: {barcodes_path}")

        codebook_path = self._resolve_codebook_path()
        if codebook_path is None:
            errors.append("No codebook configured (merlin.codebook_template).")
        elif not codebook_path.exists():
            errors.append(f"Codebook file not found: {codebook_path}")

        try:
            import plotly  # noqa: F401
        except ImportError:
            errors.append(
                "plotly is required for spatial_visualization. "
                "Install with: pip install 'merfish-pipeline[viz]'"
            )

        return errors

    def check_outputs_exist(self) -> bool:
        out = self.get_output_dir()
        return (out / "spatial_3d.html").exists()

    def run(self, dry_run: bool = False) -> StageResult:
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        cfg = self.config.spatial_visualization

        barcodes_path = self._resolve_barcodes_path()
        codebook_path = self._resolve_codebook_path()

        if barcodes_path is None or codebook_path is None:
            return StageResult(
                status="failed",
                error="Could not resolve barcodes or codebook path.",
            )

        if dry_run:
            self.logger.info(
                "[DRY RUN] Would generate 3D visualization from %s",
                barcodes_path,
            )
            return StageResult(status="skipped", metadata={"dry_run": True})

        # Check plotly availability
        try:
            import plotly  # noqa: F401
        except ImportError:
            self.logger.warning(
                "plotly not installed. Skipping spatial visualization. "
                "Install with: pip install 'merfish-pipeline[viz]'"
            )
            return StageResult(
                status="completed",
                metadata={"skipped_reason": "plotly not installed"},
            )

        # Load data
        self.logger.info("Loading barcodes from %s ...", barcodes_path)
        barcodes = read_sheet(barcodes_path)
        self.logger.info("Loaded %d barcodes.", len(barcodes))

        codebook = _load_codebook(codebook_path)

        # Merge gene symbols
        if "gene_symbol" not in barcodes.columns and "barcode_id" in barcodes.columns:
            gene_map = codebook.set_index("barcode_id")["gene_symbol"]
            barcodes["gene_symbol"] = barcodes["barcode_id"].map(gene_map)

        # Detect columns
        try:
            x_col = _detect_column(barcodes, _GLOBAL_X_CANDIDATES, "x/global_x")
            y_col = _detect_column(barcodes, _GLOBAL_Y_CANDIDATES, "y/global_y")
            z_col = _detect_column(barcodes, _Z_CANDIDATES, "z")
            fov_col = _detect_column(barcodes, _FOV_CANDIDATES, "fov")
        except ValueError as exc:
            return StageResult(status="failed", error=str(exc))

        has_cell_id = "Cell_ID" in barcodes.columns and barcodes["Cell_ID"].notna().any()

        # Downsample if configured
        if cfg.max_points is not None and len(barcodes) > cfg.max_points:
            self.logger.info(
                "Downsampling from %d to %d points (max_points=%d)",
                len(barcodes), cfg.max_points, cfg.max_points,
            )
            barcodes = barcodes.sample(n=cfg.max_points, random_state=42)

        # Build figure
        self.logger.info(
            "Building 3D figure: %d FOVs, %d barcodes, cell_id=%s",
            barcodes[fov_col].nunique(), len(barcodes), has_cell_id,
        )
        fig = _build_figure(
            barcodes, x_col, y_col, z_col, fov_col, has_cell_id, cfg,
        )

        # Write HTML
        output_dir.mkdir(parents=True, exist_ok=True)
        html_path = output_dir / "spatial_3d.html"
        fig.write_html(str(html_path), include_plotlyjs="cdn")
        self.logger.info("Wrote spatial visualization to %s", html_path)

        output_files = [str(html_path)]

        metadata = {
            "n_barcodes": len(barcodes),
            "n_fovs": int(barcodes[fov_col].nunique()),
            "has_cell_id": has_cell_id,
            "x_col": x_col,
            "y_col": y_col,
            "z_col": z_col,
        }

        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata=metadata,
        )

        self.write_run_metadata(result, start_time, parameters=metadata)
        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_barcodes_path(self) -> Path | None:
        """Resolve barcodes CSV: explicit → cell_assignment → filter_barcodes → MERlin."""
        explicit = self.config.spatial_visualization.barcodes_file
        if explicit is not None:
            return Path(explicit)

        output_dir = Path(self.config.paths.output_dir)

        ca_output = output_dir / "cell_assignment" / "barcodes_assigned.csv"
        if ca_output.exists():
            return ca_output

        fb_output = output_dir / "filter_barcodes" / "barcodes_filtered.csv"
        if fb_output.exists():
            return fb_output

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
