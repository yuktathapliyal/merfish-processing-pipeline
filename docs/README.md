# merFISH Processing Pipeline

## What is merFISH?

Multiplexed Error-Robust Fluorescence In Situ Hybridization (merFISH) is a
spatial transcriptomics technology that images hundreds of RNA species directly
inside tissue. A microscope captures thousands of images across multiple
imaging rounds, field-of-view (FOV) positions, z-slices, and fluorescent
channels. A decoding algorithm (MERlin) then converts these raw images into
spatially resolved RNA transcript counts.

## What does this pipeline do?

This pipeline handles everything **before and after** MERlin:

```
Raw Images --> [Pre-processing & QC] --> [Reregistration & Conversion] --> MERlin --> [Post-MERlin Analysis] --> AnnData
```

- **Pre-processing & QC** -- Scan raw data, build tile mosaics, detect best
  focus, and measure stage drift.
- **Reregistration & Conversion** -- Optionally correct z-drift, convert raw
  images into MERlin's expected format, and generate MERlin parameter files.
- **MERlin** -- Runs externally (not part of this pipeline). Decodes barcodes
  from the converted images.
- **Post-MERlin Analysis** -- Filter barcodes, validate against bulk RNA-seq,
  segment cells, assign barcodes to cells, generate QC reports, and export
  cell-by-gene matrices.

## Supported microscopes

| Microscope | Raw format | Convert stage |
|------------|-----------|---------------|
| **ONI** | Per-channel TIFF directories | `convert` |
| **NIKON** | Per-channel TIFF directories | `convert` |
| **ANDOR** | HDF5 `.ims` files | `ims_convert` |

## The 14 stages at a glance

| # | Stage | Phase | What it does |
|---|-------|-------|-------------|
| 1 | `index` | Pre-processing | Scan raw data, build manifest and position files |
| 2 | `stitch` | Pre-processing | Build tile mosaics for visual QC |
| 3 | `focus_qc` | Pre-processing | Find best-focus z-slice per FOV |
| 4 | `inspect_positions` | Pre-processing | Measure stage drift between rounds |
| 5 | `reregistration` | Reregistration & Conversion | Remap z-slices to uniform depth (optional) |
| 6 | `convert` | Reregistration & Conversion | Merge TIFFs for MERlin (ONI/NIKON) |
| 7 | `ims_convert` | Reregistration & Conversion | Convert IMS to merged TIFFs (ANDOR) |
| 8 | `merlin_config` | Reregistration & Conversion | Generate MERlin parameter files and launch script |
| 9 | `filter_barcodes` | Post-MERlin | Remove duplicate barcodes from reregistration |
| 10 | `correlation` | Post-MERlin | Validate barcode counts against bulk RNA-seq |
| 11 | `segmentation` | Post-MERlin | Cellpose cell segmentation on aligned images |
| 12 | `cell_assignment` | Post-MERlin | Assign decoded barcodes to segmented cells |
| 13 | `barcode_qc` | Post-MERlin | Generate QC metrics and diagnostic PDF report |
| 14 | `anndata_export` | Post-MERlin | Export cell-by-gene matrix to AnnData h5ad |

## Documentation

| Document | What it covers |
|----------|---------------|
| [Installation](installation.md) | How to install the pipeline and optional extras |
| [Configuration](configuration.md) | The config system, creating your experiment YAML, validation |
| [Pre-processing Stages](stages-preprocessing.md) | Stages 1--4: index, stitch, focus_qc, inspect_positions |
| [Reregistration & Conversion Stages](stages-reregistration.md) | Stages 5--8: reregistration, convert, ims_convert, merlin_config |
| [Post-MERlin Stages](stages-post-merlin.md) | Stages 9--14: filter_barcodes through anndata_export |
| [Workflows](workflows.md) | End-to-end examples for ONI, NIKON, and ANDOR |
| [Understanding Your Outputs](outputs-guide.md) | What each output file contains and how to interpret it |
| [CLI Reference](cli-reference.md) | All commands, flags, and examples |
| [Troubleshooting](troubleshooting.md) | Common errors and how to fix them |

## Example Data

Full worked examples with real data (raw images, pipeline outputs, MERlin
results) for all three microscope types are available on the Numbers server:

```
/projects/molonc/scratch/ythapliyal/MERFISH_EXAMPLE_FOLDER
```

These files are too large to include in the git repository.
