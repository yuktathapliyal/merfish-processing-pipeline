# Understanding Your Outputs

This guide explains what each output file contains, what to look for, and what
"good" vs "bad" results look like.

> **Full worked examples** with real data (raw images, pipeline outputs, MERlin
> results) for all three microscope types are available on the Numbers server:
>
> `/projects/molonc/scratch/ythapliyal/MERFISH_EXAMPLE_FOLDER`
>
> These files are too large to include in the git repository. Browse the example
> folder to see what a complete pipeline run looks like for ONI, NIKON, and
> ANDOR experiments.

---

## Output directory layout

After a full run, your output directory looks like this:

```
output_dir/
├── index/
│   ├── manifest.csv
│   ├── positions.standardized.csv
│   └── run_metadata.json
├── stitch/
│   └── raw/
│       └── *.TIFF
├── focus_qc/
│   ├── best_focus_slices.csv
│   ├── heatmap.png
│   └── summary.txt
├── inspect_positions/
│   ├── drift_report.csv
│   ├── drift_summary.txt
│   ├── drift_plot.png
│   └── trajectory_plot.html
├── reregistration/                (if enabled)
│   └── zmap_new_to_old.csv
├── remapped_data/                 (if reregistration ran)
│   └── {channel}/*.TIFF
├── merlin_data/                   (merged TIFFs for MERlin)
│   └── merFISH_merged_*.tiff
├── parameters/                    (MERlin parameter files)
│   ├── dataorganization/
│   ├── microscope/
│   ├── positions/
│   ├── analysis/
│   └── codebooks/
├── merlin_analysis/               (MERlin working directory)
│   ├── .merlinenv
│   ├── run_merLIN.sh
│   └── {experiment}/              (created by MERlin)
│       └── ExportBarcodes/barcodes.csv
├── merlin_config/
│   └── run_metadata.json
├── filter_barcodes/               (if enabled)
│   └── barcodes_filtered.csv
├── correlation/                   (if enabled)
│   ├── info_{name}.csv
│   ├── merged_counts/
│   └── correlation_plots.pdf
├── optimize_correlation/          (if enabled)
│   ├── correlation_trend.csv
│   ├── correlation_trend.png
│   ├── optimal_genes.csv
│   └── detailed_results.xlsx
├── segmentation/                  (if enabled)
│   ├── preprocessed/
│   └── masks/
├── cell_assignment/               (if enabled)
│   ├── barcodes_assigned.csv
│   ├── barcodes_assigned_filtered.csv
│   └── assignment_summary.csv
├── barcode_qc/                   (if enabled)
│   ├── qc_summary.csv
│   ├── per_fov_stats.csv
│   ├── per_gene_stats.csv
│   ├── per_bit_stats.csv
│   ├── per_cell_stats.csv
│   ├── qc_report.pdf
│   └── spatial_plots/
│       └── {name}_FOV_*.pdf
├── anndata_export/               (if enabled)
│   ├── {experiment}.h5ad
│   ├── cell_gene_matrix.csv
│   └── cell_metadata.csv
└── spatial_visualization/        (if enabled)
    └── spatial_3d.html
```

Every stage writes a `run_metadata.json` file in its output directory with
timestamp, status, runtime, and a list of output files.

---

## Pre-processing outputs

These outputs help you assess raw data quality before making any changes.

### `index/manifest.csv`

The master inventory of your raw data. Each row represents one image file.

| Column | Meaning |
|--------|---------|
| `round` | Imaging round number (1, 2, 3, ...) |
| `fov` | Field-of-view number |
| `z_slice` | Z-slice position in the stack |
| `channel` | Fluorescent channel name (e.g. "488nm, Raw") |
| `path` | Full file path to the raw image |

**What to check:** Count unique values in each column. Do the numbers match
your experiment design? Missing rows usually mean the pipeline couldn't find
files matching the expected naming pattern.

### `stitch/raw/*.TIFF`

Mosaic images showing all FOV tiles arranged in their correct spatial positions.

**How to view:** Open in ImageJ/FIJI. Zoom out to see the full tissue. Each
bright patch is one FOV.

**Good:** Tiles fit together seamlessly. Tissue is visible across the mosaic.
Adjacent tiles show matching features at their borders.

**Bad:** Tiles are scattered randomly, overlap incorrectly, or have dark gaps.
This usually means the position file is wrong.

### `focus_qc/heatmap.png`

A color-coded grid: rows = FOVs, columns = imaging rounds. The color indicates
the best-focus z-slice number.

**Good:** Mostly uniform colors (all FOVs focused at similar z-slices).

**Bad:** A rainbow of colors across FOVs. This means the microscope focus
drifted and you should enable reregistration.

### `inspect_positions/drift_plot.png`

Three panels showing how much the microscope stage shifted between rounds.

- **Panel 1 (line plot):** Displacement over time. Should be close to 0.
- **Panel 2 (strip plot):** Distribution of drift per round. Tight clusters = good.
- **Panel 3 (heatmap):** Displacement by FOV and round. Helps identify specific
  FOVs or rounds with problems.

**Good:** Displacement consistently < 1 pixel.

**Bad:** Large spikes or systematically increasing drift.

---

## Reregistration outputs

### `reregistration/zmap_new_to_old.csv`

Shows how z-slices were remapped for each FOV.

| Column | Meaning |
|--------|---------|
| `fov` | FOV number |
| `new_z` | New z-slice index (after remapping) |
| `old_z` | Original z-slice index (before remapping) |
| `is_duplicate` | `True` if this slice was copied from an adjacent slice to fill a gap |

**What to check:** Look at `is_duplicate`. A few `True` values are normal. If
most values are `True`, the original z-ranges were very different and
`filter_barcodes` will need to remove many duplicates after MERlin.

### `reregistration/run_metadata.json`

Contains `target_z` -- the uniform z-depth all FOVs were mapped to. This is
the minimum usable z-range across all FOVs after accounting for focus drift.

---

## MERlin parameter files

### `parameters/dataorganization/data_organization_*.csv`

This is the "decoder ring" that tells MERlin which image frame corresponds to
which bit. Each row maps a frame in the stacked TIFF to an imaging round,
channel, and color.

**What to check:** The number of rows should equal (rounds x channels x z-slices).
If this doesn't match, the data organization template may need updating.

### `merlin_analysis/run_merLIN.sh`

The launch script for MERlin. Read it before running to confirm paths are correct.

---

## Post-MERlin outputs

### `filter_barcodes/barcodes_filtered.csv`

Same format as MERlin's `barcodes.csv` but with duplicate-z barcodes removed.

**What to check:** Read `removal_summary.txt`. Normal: < 10% barcodes removed.
High removal rate (> 30%) suggests the reregistration had extensive z-duplication.

### `correlation/info_*.csv`

One row per distance threshold tested.

| Column | Meaning |
|--------|---------|
| `Distance Threshold` | How close a barcode must be to its codebook entry to be counted |
| `Pearson correlation` | Pearson r between merFISH counts and bulk RNA-seq. Higher = better. |
| `Spearman correlation` | Spearman rho (rank correlation). More robust to outliers. |
| `# detected barcodes` | How many barcodes pass this threshold |
| `# detected control barcodes` | How many blank (control) barcodes pass |

**Good:** Pearson > 0.5 at the optimal threshold. Blank count is much lower
than coding barcode count.

**Bad:** Pearson < 0.2 suggests poor decoding quality.

### `correlation/correlation_plots.pdf`

Scatter plots: each point is a gene. X-axis = bulk RNA-seq expression,
Y-axis = merFISH barcode count. Both on log2 scale.

**Good:** Points cluster along the diagonal. Blank genes (if labeled) are near
the origin.

**Bad:** Points are scattered with no correlation. Many genes at zero in one
axis but not the other.

### `optimize_correlation/correlation_trend.png`

A line plot showing the best Pearson correlation achieved at each gene group
size. The x-axis is group size, the y-axis is correlation. A red dashed line
marks your configured threshold.

**Good:** The curve rises steeply at small group sizes and plateaus around your
full panel size. This means most genes in your panel contribute positively to
the correlation.

**Bad:** The curve peaks at a small group size (e.g. 20--30 genes) then drops
as more genes are added. The genes beyond the peak are hurting your correlation
-- check `optimal_genes.csv` to see which genes perform well and which don't.

### `optimize_correlation/optimal_genes.csv`

The gene list for the best-performing group. Each row has the gene name, its
merFISH count (log-transformed), and bulk expression (log-transformed).

Use this to compare against your full panel -- genes missing from this list
may have probe issues or may not be well-expressed in this particular tissue.

### `barcode_qc/spatial_plots/*.pdf`

One PDF per FOV. Each PDF has a grid of scatter plots (one per z-slice)
showing barcode positions colored by distance to their codebook entry. The
colormap goes from dark red (close match, high confidence) to white (far
from codebook, low confidence).

**Good:** Barcodes are spread across the FOV with mostly dark colors. The
spatial distribution looks uniform.

**Bad:** A cluster of white/light points in one region suggests a localized
imaging problem (e.g. out-of-focus corner, illumination gradient). Empty
z-slices at the top or bottom may indicate the z-stack doesn't cover the
full tissue thickness.

### `spatial_visualization/spatial_3d.html`

An interactive 3D scatter plot you open in a web browser. Each point is a
decoded barcode placed at its (x, y, z) coordinates. Use the dropdown to
switch between FOVs and between gene-colored and cell-colored views. Hover
over points to see gene names and cell IDs.

**How to use:** Open the HTML file in Chrome or Firefox. Click and drag to
rotate. Scroll to zoom. Use the dropdown in the top-left to pick a FOV
and coloring mode.

### `segmentation/masks/fov_*_masks.tif`

Cell segmentation masks. Pixel values: 0 = background, 1+ = cell IDs.

**How to view:** Open in ImageJ. Apply a color lookup table (Analyze > Color >
"glasbey" or "16 colors") to see individual cells as different colors.

**Good:** Cells are cleanly separated. Most tissue area is covered. Cell sizes
look biologically reasonable.

**Bad:** Cells are merged into large blobs (under-segmentation), fragmented
into tiny pieces (over-segmentation), or most of the tissue is background
(poor detection).

### `cell_assignment/assignment_summary.csv`

Per-FOV table showing assignment rates.

| Column | Meaning |
|--------|---------|
| `fov` | FOV number |
| `n_barcodes` | Total barcodes in this FOV |
| `n_assigned` | Barcodes that landed inside a cell |
| `n_cells` | Number of unique cells found |

**What to check:** The ratio `n_assigned / n_barcodes` tells you what fraction
of barcodes were successfully assigned. Typical values: 30--70% (some barcodes
naturally fall outside cells). Consistently low assignment (< 20%) suggests
segmentation problems.

### `barcode_qc/qc_summary.csv` — SNR and error-rate columns

In addition to the basic QC columns (total_barcodes, unique_genes,
blank_barcode_pct, etc.), the summary includes these columns when per-bit
intensity data is available:

| Column | Meaning |
|--------|---------|
| `snr_contrast_median` | Median ON/OFF contrast ratio across all barcodes. Range: −1 to +1 (higher = better). Above 0.5 is good, above 0.7 is excellent. |
| `snr_contrast_mean` | Mean ON/OFF contrast ratio. |
| `snr_contrast_std` | Standard deviation of contrast. Low std = consistent quality. |
| `misid_rate` | Misidentification rate estimated from blank barcodes. Should be below 0.05 (5%). Computed even without per-bit intensities. |

See [stages-post-merlin.md — Understanding the SNR and error metrics](stages-post-merlin.md#understanding-the-snr-and-error-metrics)
for the full explanation of how these are defined and why.

### `barcode_qc/per_fov_stats.csv` — per-FOV SNR columns

| Column | Meaning |
|--------|---------|
| `snr_contrast_median` | Median contrast ratio for barcodes in this FOV |
| `misid_rate` | Per-FOV misidentification rate |

Compare these across FOVs to spot spatial patterns in signal quality. A FOV
with much lower SNR or higher misidentification rate than others may have had
an imaging problem.

### `barcode_qc/per_bit_stats.csv`

One row per bit (imaging round). Helps identify weak readout sequences.

| Column | Meaning |
|--------|---------|
| `bit_index` | Bit position (0-indexed, matches `intensity_0`, `intensity_1`, ...) |
| `bit_name` | Readout sequence name from the codebook (e.g. RS0015) |
| `median_on` | Median intensity for barcodes where this bit should be ON |
| `median_off` | Median intensity for barcodes where this bit should be OFF |
| `contrast` | `(median_on − median_off) / (median_on + median_off)` |
| `n_on` | Number of barcodes contributing to the ON median |
| `n_off` | Number of barcodes contributing to the OFF median |

**Good:** All bits have similar contrast values (e.g. 0.85 -- 0.92).

**Bad:** One or two bits with much lower contrast than the rest (e.g. 0.4
when others are 0.8+). That readout sequence likely had a problem -- check the
raw images for that imaging round.

---

## The QC report (`barcode_qc/qc_report.pdf`)

This is the most important output for quickly assessing experiment quality. It
contains up to 8 diagnostic panels (6 base panels + 2 SNR panels when per-bit
intensity data is available):

### Panel 1: Barcode abundance (top-left)

Genes ranked by barcode count (log scale). Coding genes in blue, blank control
genes in red.

**Good:** Smooth descending curve for coding genes. Blanks are at the
bottom-right with much lower counts.

**Bad:** Flat distribution (all genes have similar counts) suggests poor
decoding. Blanks mixed in with coding genes suggests high noise.

### Panel 2: Intensity distribution (top-right)

Histogram of log10(mean_intensity) for all barcodes.

**Good:** Single clear peak. The peak position depends on your imaging
conditions.

**Bad:** Bimodal (two peaks) may indicate a mix of real barcodes and noise.
Very broad distribution suggests inconsistent imaging.

### Panel 3: Barcodes per FOV (middle-left)

Bar chart showing barcode count per FOV.

**Good:** Roughly uniform across FOVs. Some variation is normal depending on
tissue density.

**Bad:** Some FOVs with dramatically fewer barcodes may have imaging problems
(out of focus, low signal, etc.).

### Panel 4: Distance threshold curve (middle-right)

Pearson and Spearman correlation vs distance threshold. Only appears if
MERlin's `PlotPerformance/` directory exists.

**Good:** Curves peak at an intermediate threshold (typically 0.3--0.5).
The vertical red line marks the optimal threshold.

**Bad:** Flat curves near zero at all thresholds.

### Panel 5: Barcodes per cell (bottom-left)

Histogram showing how many barcodes each cell received. Only appears if
`cell_assignment` was run.

**Good:** Right-skewed distribution with most cells having 20--200 barcodes.
Median (red line) in the tens to hundreds range.

**Bad:** Most cells having 0--5 barcodes suggests very sparse decoding or
poor segmentation.

### Panel 6: Top genes (bottom-right)

Horizontal bar chart of the N most abundant genes.

**What to check:** Do the top genes make biological sense for your tissue type?
Housekeeping genes often appear near the top.

### Panel 7: SNR contrast distribution

Histogram of the per-barcode ON/OFF contrast ratio. Only appears if the
barcodes CSV has per-bit intensity columns (`intensity_0`..`intensity_N`) and
a codebook is configured. The red dashed line marks the median. Colored
vertical lines mark the quality thresholds: orange (0.3 = poor/marginal
boundary), gold (0.5 = marginal/good), green (0.7 = good/excellent).

**Good:** A narrow peak well above 0.5, ideally above 0.7. The further right
the peak, the better the ON/OFF separation.

**Bad:** A broad distribution centered below 0.5, or a bimodal shape with
a second peak near 0. The low-contrast barcodes likely have decoding errors.

### Panel 8: Per-bit ON/OFF contrast

Bar chart showing the contrast ratio for each bit (imaging round). Bars are
color-coded: green (>= 0.7), gold (0.5 -- 0.7), orange (0.3 -- 0.5), red
(< 0.3). A gray dashed line marks 0.5.

**Good:** All bars are similar height and green/gold.

**Bad:** One or two bars are orange/red while the rest are green. That readout
sequence had problems -- check the raw images for that imaging round. Possible
causes: weak hybridisation, poor focus during that round, a bad readout probe.

---

## Final export

### `anndata_export/{experiment}.h5ad`

The AnnData object -- your entry point to single-cell analysis in Python.

```python
import anndata
adata = anndata.read_h5ad("path/to/experiment.h5ad")

print(adata.shape)              # (n_cells, n_genes)
print(adata.obs.head())         # cell metadata: fov, n_barcodes, n_genes
print(adata.var.head())         # gene metadata: is_blank
print(adata.obsm['spatial'])    # spatial coordinates (global_x, global_y)
```

### `anndata_export/cell_gene_matrix.csv`

Plain CSV: rows = cells, columns = genes. Works without anndata installed.
Open in Excel for quick inspection or load with `pd.read_csv()`.

### `anndata_export/cell_metadata.csv`

Per-cell metadata. Use to filter low-quality cells:

```python
import pandas as pd
meta = pd.read_csv("cell_metadata.csv", index_col=0)
# Filter cells with fewer than 10 barcodes
good_cells = meta[meta['n_barcodes'] >= 10].index
```

---

Back to: [Documentation Index](README.md)
