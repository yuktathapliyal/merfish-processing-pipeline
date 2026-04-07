# Post-MERlin Analysis (Stages 9--16)

These stages run **after MERlin has finished decoding barcodes**. They operate
on the `ExportBarcodes/barcodes.csv` file that MERlin produces -- a table with
one row per detected barcode, containing its spatial coordinates, intensity,
distance from the nearest codebook entry, FOV, and barcode ID.

The post-MERlin stages filter, validate, and enrich this barcode data:
filtering out reregistration artifacts, checking decoding quality against bulk
RNA-seq, identifying which barcodes belong to which cells, generating QC
reports, and finally exporting a cell-by-gene count matrix for downstream
single-cell analysis.

**Auto-detection:** Most stages in this group automatically find their inputs
from upstream stage outputs. For example, `cell_assignment` looks for filtered
barcodes from `filter_barcodes` first, then falls back to MERlin's raw output.
You rarely need to set explicit file paths.

---

## `filter_barcodes`

**What it does:** If you ran `reregistration`, some z-slices were duplicated
(copied from adjacent slices to fill gaps). MERlin doesn't know about this
duplication, so it may decode the same RNA molecule twice -- once from the
original z-slice and once from its duplicate. This stage identifies and removes
those duplicate barcodes.

**When to use:** Only if you enabled `reregistration`. If you didn't use
reregistration, skip this stage.

**Output:**

| File | What it contains |
|------|-----------------|
| `filter_barcodes/barcodes_filtered.csv` | Same columns as MERlin's `barcodes.csv` but with duplicate-z barcodes removed. This is the "clean" barcode table that all downstream stages should use. |
| `filter_barcodes/removal_summary.txt` | How many barcodes were removed and why. A summary showing total barcodes before/after filtering and the number of duplicates found. |

**What to check:** Read `removal_summary.txt`. A small percentage of removed
barcodes (< 10%) is normal. If a large fraction was removed, check your
reregistration z-map for excessive duplication.

**Config:**

```yaml
filter_barcodes:
  enabled: true
  mode: "any"           # "any" = remove if ANY round has a duplicate at this z
                        # "all" = only remove if ALL rounds have duplicates (stricter)
  barcodes_file: null   # auto-detected from MERlin output, or set a path explicitly
```

**Requires:** `reregistration` completed + MERlin completed externally.

---

## `correlation`

**What it does:** Validates the quality of barcode decoding by comparing the
decoded barcode counts against an independent bulk RNA-seq expression dataset.
If the merFISH experiment worked well, the gene expression levels measured by
merFISH should correlate with the bulk RNA-seq measurements.

For each distance threshold (barcodes closer to their codebook entry are more
likely to be real), this stage counts barcodes per gene, merges with the bulk
data, and computes Pearson and Spearman correlation coefficients.

**Output:**

| File | What it contains |
|------|-----------------|
| `correlation/info_{name}.csv` | Summary table with one row per distance threshold: threshold value, Pearson correlation, Spearman correlation, number of detected barcodes, number of blank (control) barcodes, and number of genes detected. |
| `correlation/merged_counts/{name}_{threshold}.csv` | Per-gene table at each threshold: gene name, merFISH barcode count, bulk RNA-seq expression value. Useful for finding outlier genes. |
| `correlation/correlation_plots.pdf` | Scatter plots showing merFISH vs bulk expression for each threshold. Each point is a gene. Points should cluster along the diagonal if decoding worked well. |

**What to check:** Look at `info_{name}.csv`. **Good experiments** typically
show Pearson correlation > 0.5 at the optimal distance threshold. The optimal
threshold is the one with the highest Pearson value. Low correlations (< 0.2)
suggest problems with the experiment or imaging quality.

In the scatter plots, genes that fall far from the diagonal may have probe
issues or be poorly expressed in your tissue.

**Config:**

```yaml
correlation:
  enabled: true
  bulk_file: "/path/to/bulk_expression.csv"     # required: must have gene_symbol and expression columns
  distance_thresholds: [0.5167, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25]
```

**Requires:** MERlin completed externally + codebook configured
(`merlin.codebook_template`).

---

## `optimize_correlation`

**What it does:** Takes the gene-level correlation data from the `correlation`
stage and searches for the subset of genes that gives the best possible
correlation with your bulk RNA-seq reference. It uses a technique called
simulated annealing -- a randomized optimization algorithm that tries swapping
genes in and out of a group to maximize the Pearson correlation.

This is useful for identifying which genes in your panel are performing well
and which are dragging down the overall correlation. It also helps determine
the right panel size -- the stage tests group sizes from small (5 genes) up to
large (120 genes) and shows how correlation changes with group size.

**When to use:** After running `correlation`, when you want to understand which
genes are contributing to or hurting your experiment's correlation with bulk
data. Particularly useful for panel design and troubleshooting low-performing
experiments.

**Output:**

| File | What it contains |
|------|-----------------|
| `optimize_correlation/correlation_trend.csv` | One row per group size tested: the size and the best Pearson correlation achieved at that size. |
| `optimize_correlation/correlation_trend.png` | Line plot showing how correlation changes with group size. The red dashed line marks your correlation threshold. |
| `optimize_correlation/optimal_genes.csv` | The gene list for the best-performing group: gene symbol, log-transformed merFISH counts, log-transformed bulk expression, and the group's correlation value. |
| `optimize_correlation/detailed_results.xlsx` | Excel workbook with a summary sheet and one sheet per group size, listing the genes in each optimized group. |

**What to check:** Look at `correlation_trend.png` first. Correlation typically
increases with group size up to a point, then plateaus or drops as weaker genes
get included. The peak tells you the effective panel size. Check
`optimal_genes.csv` to see which genes made it into the best group -- genes
that never appear in any optimized group may have probe issues.

**Config:**

```yaml
optimize_correlation:
  enabled: true
  correlation_threshold: 0.45   # minimum Pearson r to accept a group
  n_attempts: 5                 # independent optimization runs per size
  # size_range_start: 5         # smallest group size to test
  # size_range_end: 120         # largest group size to test
  # size_range_step: 5          # step between sizes
  # max_iterations: 2000        # SA iterations per attempt
  # cooling_rate: 0.995         # temperature decay rate
  # distance_threshold: null    # which correlation threshold's data to use (auto-selects best)
  # random_seed: null           # set for reproducible results
```

**Requires:** `correlation` stage completed (reads its `merged_counts/` output).

---

## `segmentation`

**What it does:** Identifies individual cells in the tissue using
[Cellpose](https://github.com/MouseLand/cellpose), a deep-learning cell
segmentation model. This is needed so that decoded RNA barcodes can be assigned
to specific cells (rather than just floating in space).

The stage preprocesses MERlin's aligned images to create a two-channel input
(nuclei + cytoplasm), then runs Cellpose to produce segmentation masks --
images where each pixel is labeled with the ID of the cell it belongs to (or 0
for background).

**How it works:**

1. **Preprocessing** (per FOV): Loads an aligned multi-frame TIFF stack from
   MERlin. Extracts a **nuclei channel** (one bit, set by `nuclei_bit`) and
   builds a **cytoplasm channel** by summing all remaining bits. Applies a
   median filter and min-max normalizes both channels to [0, 1]. Output is a
   4-D volume `(Z, 2, Y, X)` where channel 0 = cytoplasm, channel 1 = nuclei.

2. **Segmentation** (per FOV): Runs Cellpose in **2D+stitch mode** -- each
   z-slice is segmented independently as a 2D image, then masks are stitched
   across z-slices using an IoU overlap threshold (`stitch_threshold`). If a
   cell in slice z=5 overlaps sufficiently with a cell in z=6, they receive
   the same cell ID.

**Output:**

| File | What it contains |
|------|-----------------|
| `segmentation/preprocessed/fov_*_preprocessed.tif` | Preprocessed 4-D volumes with nuclei and cytoplasm channels separated. Useful for debugging if segmentation results look wrong. |
| `segmentation/masks/fov_*_masks.tif` | Cell segmentation masks. Each pixel value is a cell ID (1, 2, 3, ...) or 0 for background. Open in ImageJ to see colored cell boundaries. |

**What to check:** Open a mask TIFF in ImageJ (use LUT > "glasbey" for
colors). Each cell should be a distinct colored blob. Check that cells are
properly separated (not merged) and that most of the tissue is segmented (not
too much background). Small fragments or very large merged blobs suggest the
`diameter` parameter needs adjustment.

**Config (3D mode, default):**

```yaml
segmentation:
  enabled: true
  aligned_images_dir: "/path/to/FiducialCorrelationWarp/images"  # required
  mode: "3d"                # segment all z-slices with 2D+stitch (default)
  nuclei_bit: 17            # 0-indexed bit for nuclei channel
  total_bits: 18            # total bits in aligned image stacks
  exclude_bits: []          # bits to exclude from cytoplasm (e.g. [0] for fiducial)
  median_kernel: 3          # median filter kernel size (odd integer)
  model_type: "cpsam"       # Cellpose model ("cyto2", "cpsam", etc.)
  diameter: null             # cell diameter in pixels (null = auto-detect)
  batch_size: 8             # GPU batch size
  stitch_threshold: 0.5     # IoU threshold for stitching masks across z
```

**Config (2D mode -- single z-slice):**

```yaml
segmentation:
  enabled: true
  aligned_images_dir: "/path/to/FiducialCorrelationWarp/images"
  mode: "2d"                # segment only one z-slice
  reference_z_slice: 5      # which z-slice (1-indexed by default)
  z_indexing: 1              # 1 = 1-indexed, 0 = 0-indexed
```

**Output format:**
- **3D mode:** Mask shape is `(Z, Y, X)` -- one label image per z-slice with
  consistent cell IDs across z.
- **2D mode:** Mask shape is `(Y, X)` -- single label image. Useful when only
  one z-slice has good enough quality for segmentation.

**Requires:** MERlin completed externally. Install with
`pip install -e ".[segmentation]"`. GPU recommended but CPU works too (slower).

---

## `cell_assignment`

**What it does:** Takes the decoded barcodes (with their pixel coordinates) and
the segmentation masks, and determines which cell each barcode belongs to. For
each barcode, it looks up the pixel position in the FOV's mask to get the cell
label. Barcodes that land on background (label 0) are marked as unassigned.

Optionally filters "border cells" -- cells with barcodes near the edge of a
FOV, which may be incomplete or duplicated across adjacent FOVs.

**Output:**

| File | What it contains |
|------|-----------------|
| `cell_assignment/barcodes_assigned.csv` | The full barcodes table (same columns as MERlin output) with a new `Cell_ID` column. Each barcode now has a cell assignment like `Cell3_42` (FOV 3, cell label 42) or `None` (background). |
| `cell_assignment/barcodes_assigned_filtered.csv` | Same table but with border cells removed. Only created when `crop_margin > 0`. |
| `cell_assignment/assignment_summary.csv` | Per-FOV summary: how many barcodes total, how many were assigned to cells, how many unique cells were found. Useful for spotting FOVs with low assignment rates. |

**What to check:** Look at `assignment_summary.csv`. The ratio of assigned vs
total barcodes tells you how well the segmentation captured the tissue. A high
unassigned rate (> 50%) may indicate poor segmentation, misaligned masks, or
lots of barcodes outside cells.

**Config:**

```yaml
cell_assignment:
  enabled: true
  # barcodes_file: null      # auto-detected from filter_barcodes or MERlin output
  # masks_dir: null           # auto-detected from {output_dir}/segmentation/masks/
  crop_margin: 0              # pixels from FOV edge for border filtering (0 = off)
```

**Auto-detection:**
- **Barcodes:** Uses `filter_barcodes/barcodes_filtered.csv` if it exists
  (preferred after reregistration), otherwise MERlin's `ExportBarcodes/barcodes.csv`.
- **Masks:** Uses `{output_dir}/segmentation/masks/` and matches FOV numbers
  from mask filenames.

**Requires:** `segmentation` stage completed + MERlin completed externally.

---

## `barcode_qc`

**What it does:** Generates a comprehensive quality control report from the
decoded barcodes. Computes summary statistics, per-gene/FOV/cell metrics, and
creates a multi-panel PDF with diagnostic plots. This automates the manual QC
analysis you'd otherwise do in a Jupyter notebook.

If the barcodes have a `Cell_ID` column (from `cell_assignment`), it also
includes per-cell statistics. If MERlin's `PlotPerformance/` directory exists,
it reads the distance-threshold correlation data from there.

**Output:**

| File | What it contains |
|------|-----------------|
| `barcode_qc/qc_summary.csv` | A single-row table with all key metrics: total barcodes, unique genes, blank barcode rate, barcodes per FOV stats, barcodes per cell stats, optimal distance threshold, and correlation values. One row = one experiment summary. |
| `barcode_qc/per_fov_stats.csv` | Barcode count and mean intensity per FOV. Useful for identifying FOVs with unusually low or high barcode counts. |
| `barcode_qc/per_gene_stats.csv` | Barcode count per gene, sorted by abundance, with an `is_blank` flag. Blank genes (controls) should have much fewer barcodes than coding genes. |
| `barcode_qc/per_cell_stats.csv` | Barcodes per cell and genes per cell. Only created if `Cell_ID` is present in the barcodes. |
| `barcode_qc/qc_report.pdf` | A 6-panel diagnostic plot. See the [outputs guide](outputs-guide.md) for how to read each panel. |
| `barcode_qc/spatial_plots/{name}_FOV_NNN.pdf` | One PDF per FOV showing barcode positions colored by distance to codebook. Each subplot is a different z-slice. Helps identify spatial patterns in decoding quality. |

**What to check:** Start with `qc_report.pdf` -- it gives you a visual
overview of experiment quality in one page. Then check `qc_summary.csv` for
the numbers. Key metrics to look at:
- `blank_barcode_pct` -- should be low (< 5% for a good experiment)
- `barcodes_per_fov_cv` -- coefficient of variation; low = uniform across FOVs
- `barcodes_per_cell_median` -- typical range: 20--200 depending on tissue

The per-FOV spatial plots in `spatial_plots/` show where barcodes are landing
within each FOV. Red spots indicate barcodes far from their codebook entry
(low confidence). Look for spatial patterns -- a corner of the FOV that's
consistently red may indicate an optical issue.

**Config:**

```yaml
barcode_qc:
  enabled: true
  top_n_genes: 20              # how many top genes to show in the report
  spatial_plots_enabled: true  # generate per-FOV scatter plots (one PDF each)
  spatial_plots_columns: 3     # columns in the z-slice subplot grid
  # barcodes_file: null        # auto-detected
```

**Auto-detection (barcodes):**
1. `cell_assignment/barcodes_assigned.csv` (preferred -- includes Cell_ID)
2. `filter_barcodes/barcodes_filtered.csv`
3. MERlin `ExportBarcodes/barcodes.csv`

**Requires:** MERlin completed externally + codebook configured.

---

## `anndata_export`

**What it does:** This is the final pipeline stage. It takes the cell-assigned
barcodes and builds a **cell-by-gene count matrix** -- a table where rows are
cells and columns are genes, and each value is how many times that gene's
barcode was detected in that cell. This is the standard input format for
single-cell analysis tools like [scanpy](https://scanpy.readthedocs.io/) and
[squidpy](https://squidpy.readthedocs.io/).

The matrix is exported as an [AnnData](https://anndata.readthedocs.io/) h5ad
file (the standard format for single-cell data in Python) and as a plain CSV
fallback.

**Output:**

| File | What it contains |
|------|-----------------|
| `anndata_export/{experiment}.h5ad` | An AnnData object containing: the count matrix (`adata.X`), cell metadata like FOV and barcode count (`adata.obs`), gene metadata like is_blank flag (`adata.var`), and spatial coordinates (`adata.obsm['spatial']`). Load with `anndata.read_h5ad()`. |
| `anndata_export/cell_gene_matrix.csv` | Plain CSV fallback: rows = cells, columns = genes, values = barcode counts. Works without anndata installed. Open in Excel or pandas. |
| `anndata_export/cell_metadata.csv` | Per-cell metadata: Cell_ID, FOV, number of barcodes, number of unique genes. Use this to filter low-quality cells. |

**What to check:**
- Load the h5ad in Python: `import anndata; adata = anndata.read_h5ad("path.h5ad")`
- Check `adata.shape` -- should be `(n_cells, n_genes)`
- Check `adata.obs` for cell metadata, `adata.obsm['spatial']` for coordinates
- In `cell_metadata.csv`, cells with very few barcodes (< 5) may be noise

**Config:**

```yaml
anndata_export:
  enabled: true
  min_barcodes_per_cell: 0    # filter cells with fewer barcodes (0 = keep all)
  exclude_blanks: true         # remove blank (control) genes from the count matrix
  # barcodes_file: null        # auto-detected from cell_assignment
```

**Install:** The h5ad export requires the `anndata` package:
`pip install -e ".[export]"`. The CSV fallback always works without it.

**Requires:** `cell_assignment` stage completed + codebook configured.

---

## `spatial_visualization`

**What it does:** Creates an interactive 3D scatter plot of all decoded
barcodes, viewable in any web browser. Each barcode is placed at its spatial
coordinates (x, y, z), and you can rotate, zoom, and pan to explore the
tissue in three dimensions. A dropdown lets you switch between FOVs and between
two coloring modes:

- **Gene-colored:** Each barcode is colored by which gene it was decoded as.
  Useful for seeing where specific genes are expressed in the tissue.
- **Cell-colored:** Each barcode is colored by which cell it was assigned to
  (from `cell_assignment`). Unassigned barcodes are shown in dim gray. Only
  available when Cell_ID is present.

**When to use:** After all other analysis is done. This is a visualization
stage -- it doesn't produce data for downstream analysis, but it's a powerful
way to explore the spatial distribution of gene expression in your tissue.

**Output:**

| File | What it contains |
|------|-----------------|
| `spatial_visualization/spatial_3d.html` | Interactive plotly HTML file. Open in a browser to explore. Use the dropdown in the top-left to switch between FOVs and coloring modes. |

**What to check:** Rotate the 3D view to confirm barcodes are distributed
across all z-slices (not concentrated in just one). In gene-colored mode, look
for spatial clustering of specific genes -- this is expected for genes with
known spatial expression patterns. In cell-colored mode, check that cells form
coherent spatial clusters rather than a random salt-and-pepper pattern.

**Config:**

```yaml
spatial_visualization:
  enabled: true
  # marker_size: 2            # size of scatter points
  # max_points: null           # downsample to this many barcodes (null = all)
  # barcodes_file: null        # auto-detected
```

**Install:** Requires `plotly`: `pip install -e ".[viz]"`. The stage skips
gracefully if plotly is not installed.

**Requires:** MERlin completed externally + codebook configured. Best results
when `cell_assignment` has been run (enables cell-colored view).

---

Next: [Workflows](workflows.md) | [Understanding Your Outputs](outputs-guide.md)
