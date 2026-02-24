# merFISH Processing Pipeline -- User Manual

## Installation

```bash
# Clone the repository
git clone <repo-url>
cd merfish-processing-pipeline

# Install in editable mode (core pipeline, no segmentation)
pip install -e .
```

### Segmentation (cellpose) setup

**Local machine with GPU (Windows/Linux):**

```bash
# Install PyTorch with CUDA (check your CUDA version with nvidia-smi)
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121

# Install cellpose and pipeline segmentation extra
pip install -e ".[segmentation]"

# Verify GPU is available
python -c "import torch; print(torch.cuda.is_available())"  # should print True
```

**Server (CPU only):**

```bash
# Install CPU-only PyTorch (much smaller download)
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install -e ".[segmentation]"
```

**If pip fails to build fastremap** (common on servers with old GCC):

```bash
conda install -c conda-forge cellpose
```

### Optional extras

```bash
# Interactive trajectory plots (plotly)
pip install -e ".[viz]"

# Everything (segmentation + viz)
pip install -e ".[all]"
```

Verify the install:

```bash
merfish-pipe --help
```

---

## Quick Start

### 1. Generate an experiment config

```bash
merfish-pipe config init --microscope oni -o my_experiment.yaml
```

This creates an annotated YAML template. Replace `oni` with `nikon` or `andor` to
match your microscope.

### 2. Edit the required fields

Open `my_experiment.yaml` and fill in:

```yaml
experiment:
  name: "XP17596"              # your experiment ID
  microscope: "oni"             # oni, nikon, or andor

paths:
  raw_data_dir: "/data/raw/XP17596"      # where your raw images are
  output_dir: "/data/output/XP17596"     # where results go

raw_data:
  bead_channel_folder: "488nm, Raw"      # fiducial channel subfolder name

pipeline:
  stages: [index, stitch, focus_qc, inspect_positions]
```

### 3. Validate your config

```bash
merfish-pipe config validate my_experiment.yaml
```

If valid, prints "Configuration is valid." If not, shows what's wrong.

### 4. Run the pipeline

```bash
merfish-pipe run my_experiment.yaml
```

### 5. Check status

```bash
merfish-pipe status my_experiment.yaml
```

Shows which stages have completed, failed, or are pending.

---

## Configuration

The pipeline uses three config layers that merge together:

| Layer | What it controls | Location |
|-------|-----------------|----------|
| Microscope | Pixel size, image dimensions, file patterns, orientation | `configs/microscopes/{oni,nikon,andor}.yaml` |
| Experiment | Your data paths, stage settings, codebook | Your YAML file (e.g. `my_experiment.yaml`) |
| Execution | Local vs SLURM, worker count | `configs/profiles/{local,slurm}.yaml` |

**How merging works:** The microscope config provides defaults. Your experiment
config overrides what it sets. The execution profile adds runtime settings. You
only need to write the experiment config -- the other two are pre-built.

### Derived paths

These directories are auto-created under `output_dir` unless you override them:

```
output_dir/
  remapped_data/      # reregistered images (if reregistration enabled)
  merlin_analysis/    # MERlin output location
  parameters/         # MERlin parameter files
  logs/               # pipeline logs
```

### View the fully resolved config

```bash
merfish-pipe config show my_experiment.yaml
merfish-pipe config show my_experiment.yaml --json
```

---

## CLI Reference

### `merfish-pipe run`

```bash
merfish-pipe run <experiment.yaml> [OPTIONS]
```

| Flag | What it does |
|------|-------------|
| `--stage <name>` | Run only this one stage |
| `--from-stage <name>` | Run from this stage onward (skips earlier stages) |
| `--dry-run` | Show what would run without executing |
| `--force` | Re-run stages even if outputs exist |
| `--profile slurm` | Submit as SLURM job instead of running locally |
| `--workers <N>` | Override number of parallel threads |
| `-v` | Verbose logging (DEBUG level) |

Examples:

```bash
# Run all stages listed in your config
merfish-pipe run my_experiment.yaml

# Run only the index stage
merfish-pipe run my_experiment.yaml --stage index

# Resume from stitch onward
merfish-pipe run my_experiment.yaml --from-stage stitch

# Preview without running
merfish-pipe run my_experiment.yaml --dry-run

# Force re-run stitch (ignores existing output)
merfish-pipe run my_experiment.yaml --stage stitch --force

# Run on SLURM cluster
merfish-pipe run my_experiment.yaml --profile slurm
```

### `merfish-pipe config`

```bash
# Generate a new experiment config template
merfish-pipe config init --microscope oni -o my_experiment.yaml

# Validate a config file
merfish-pipe config validate my_experiment.yaml

# Show the fully merged config
merfish-pipe config show my_experiment.yaml
```

### `merfish-pipe status`

```bash
merfish-pipe status my_experiment.yaml
```

Reads `run_metadata.json` from each stage output directory and reports
completion status.

---

## Pipeline Stages

Stages run in the order listed below. Add the ones you need to your
`pipeline.stages` list.

### 1. `index`

Scans your raw data directory and builds two reference files used by all
downstream stages.

**Output:**

| File | Contents |
|------|----------|
| `index/manifest.csv` | One row per raw image: round, fov, z-slice, channel, path |
| `index/positions.standardized.csv` | Stage positions in a unified format across all microscope types |

**Config:** No special options. Just needs `paths.raw_data_dir` to point at your
raw images.

---

### 2. `stitch`

Builds tile mosaics from bead-channel (fiducial) images. Useful for visually
inspecting tissue coverage and tile overlap.

**Output:**

| File | Contents |
|------|----------|
| `stitch/raw/*.TIFF` | Stitched mosaics (one per imaging round or z-slice) |

If reregistration was run first, outputs go to `stitch/reregistered/` instead,
so raw stitches are preserved for comparison.

**Config:**

```yaml
stitch:
  group_by: "ir"        # "ir" = one mosaic per round, "z" = one per z-slice
  ir_range: [1, 9]      # optional: restrict to these rounds
  z_range: [1, 15]      # optional: restrict to these z-slices
```

**Requires:** `index` stage completed.

---

### 3. `focus_qc`

Finds the best-focus z-slice for each FOV in each imaging round using
Laplacian-of-Gaussian variance on the bead channel.

**Output:**

| File | Contents |
|------|----------|
| `focus_qc/best_focus_slices.csv` | Best z-slice per FOV per round |
| `focus_qc/heatmap.png` | Color-coded z-variation across FOVs and rounds |
| `focus_qc/summary.txt` | Aggregate focus statistics |

**Config:**

```yaml
focus_qc:
  sigma: 1.0    # Gaussian blur sigma
  ksize: 3      # kernel size
```

---

### 4. `inspect_positions`

Measures microscope stage drift between imaging rounds. Helps you decide if
reregistration is needed.

**Output:**

| File | Contents |
|------|----------|
| `inspect_positions/drift_report.csv` | Per-FOV, per-round drift (delta_x, delta_y, displacement) |
| `inspect_positions/drift_summary.txt` | Aggregate drift statistics |
| `inspect_positions/drift_plot.png` | Three-panel plot: displacement trend, strip plot by round, FOV-round heatmap |
| `inspect_positions/trajectory_plot.html` | Interactive 3D stage trajectory (needs `plotly`) |

**Config:**

```yaml
inspect_positions:
  log_file: null                # auto-detected from raw data
  rounds_to_check: null         # check all rounds by default
  trajectory_z_slices: null     # default: first 3 z-slices
```

**Requires:** `index` stage completed, **or** `ims_convert` stage completed (ANDOR
workflow -- reads per-round `stagePos_Round#N.csv` files from `merlin_data_dir`).

---

### 5. `reregistration`

Remaps z-slices so every FOV has the same depth. Needed when different FOVs
were imaged at different z-ranges due to focus drift.

Uses the best-focus data from `focus_qc` to compute a uniform target depth,
then copies and renames raw image files into `remapped_data/`.

**Output:**

| File | Contents |
|------|----------|
| `reregistration/zmap_new_to_old.csv` | Mapping of new z-indices to old, with duplicate flags |
| `remapped_data/{channel}/*.TIFF` | Remapped image files |

**Config:**

```yaml
reregistration:
  enabled: true         # must be explicitly enabled
  total_z: null         # auto-detected from data
  target_z: null        # auto-computed from best-focus
```

**Requires:** `focus_qc` stage completed.

---

### 6. `convert` (ONI and NIKON only)

Merges per-channel single-plane TIFFs into multi-frame stacked TIFFs that
MERlin expects. Each output file contains all channels and z-slices for one
(round, FOV) combination.

Automatically uses remapped data if reregistration was run.

**Output:**

| File | Contents |
|------|----------|
| `convert/merFISH_merged_{round}_{fov}.tiff` | Stacked multi-frame TIFF per round and FOV |

**Config:** No special options needed.

**Requires:** `index` stage completed. Use `ims_convert` instead for Andor data.

---

### 7. `ims_convert` (Andor only)

Converts Andor IMS (HDF5) files to the same merged TIFF format. Also extracts
stage positions from IMS metadata and writes per-round position CSVs (used
downstream by `inspect_positions` and `merlin_config`).

**Output:**

| File | Contents |
|------|----------|
| `merlin_data/merFISH_merged_{round}_{fov}.tiff` | Merged stacked TIFFs (in `merlin_data_dir`) |
| `merlin_data/stagePos_Round#{round}.csv` | Stage position CSVs (in `merlin_data_dir`) |

**Config:**

```yaml
raw_data:
  andor:
    channel_order: [0, 2, 1]    # channel reordering for your Andor setup
```

**Requires:** Raw IMS files organised in round folders (e.g. `1st round/`, `R1/`).

---

### 8. `merlin_config`

Generates all parameter files MERlin needs and a shell script to launch MERlin.
Does NOT run MERlin itself. Files are written to two directories:

**Parameter files** (saved to `parameters/`, organized by type):

| File | Contents |
|------|----------|
| `parameters/dataorganization/data_organization_{name}.csv` | Frame-to-bit mapping (expanded from template) |
| `parameters/microscope/microscope_{name}.json` | Microscope parameters for MERlin |
| `parameters/positions/positions_{name}.csv` | FOV positions for MERlin |
| `parameters/analysis/analysis_{name}.json` | Analysis task parameters |
| `parameters/codebooks/{codebook}.csv` | Barcode codebook |

**Launch files** (saved to `merlin_analysis/`):

| File | Contents |
|------|----------|
| `merlin_analysis/.merlinenv` | Sets DATA_HOME, ANALYSIS_HOME, PARAMETERS_HOME |
| `merlin_analysis/run_merLIN.sh` | Shell script to invoke MERlin |

After this stage completes, run MERlin externally:

```bash
source merlin_analysis/.merlinenv
bash merlin_analysis/run_merLIN.sh
```

**Config:**

```yaml
merlin:
  codebook_template: "/path/to/codebook.csv"   # required
  analysis_template: null                        # uses microscope default
  microscope_template: null                      # uses microscope default
  cores: 100
```

**Requires:** `convert` or `ims_convert` stage completed.

---

### 9. `filter_barcodes`

Removes barcodes from duplicated z-slices introduced by reregistration. Only
needed if you ran `reregistration`.

**Output:**

| File | Contents |
|------|----------|
| `filter_barcodes/barcodes_filtered.csv` | Cleaned barcodes without duplicates |
| `filter_barcodes/removal_summary.txt` | How many barcodes were removed and why |

Auto-detects the barcodes file at
`{output_dir}/merlin_analysis/{experiment_name}/ExportBarcodes/barcodes.csv`.

**Config:**

```yaml
filter_barcodes:
  enabled: true
  mode: "any"           # "any" = remove if ANY round has duplicate, "all" = only if ALL do
  barcodes_file: null   # auto-detected, or set an explicit path
```

**Requires:** `reregistration` completed + MERlin completed externally.

---

### 10. `correlation`

Compares decoded barcode counts against bulk RNA-seq expression data. Produces
Pearson and Spearman correlations at multiple distance thresholds.

**Output:**

| File | Contents |
|------|----------|
| `correlation/info_{name}.csv` | Summary: threshold, Pearson, Spearman, barcode counts, blank counts |
| `correlation/merged_counts/{name}_{threshold}.csv` | Per-gene counts merged with bulk expression at each threshold |
| `correlation/correlation_plots.pdf` | Scatter plots with gene labels at each threshold |

**Config:**

```yaml
correlation:
  enabled: true
  bulk_file: "/path/to/bulk_expression.csv"     # required: gene_symbol + TPM/FPKM column
  distance_thresholds: [0.5167, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25]
```

**Requires:** MERlin completed externally + codebook configured.

---

### 11. `segmentation`

Runs [Cellpose](https://github.com/MouseLand/cellpose) cell segmentation on
MERlin's aligned images. Identifies individual cells so that decoded RNA
transcripts can later be assigned to specific cells.

**How it works:**

1. **Preprocessing** (per FOV): Loads an aligned multi-frame TIFF stack from
   MERlin. Extracts a **nuclei channel** (one bit, set by `nuclei_bit`) and
   builds a **cytoplasm channel** by summing all remaining bits. Applies a
   median filter and min-max normalizes both channels to [0, 1]. Output is a
   4-D volume `(Z, 2, Y, X)` where channel 0 = cytoplasm, channel 1 = nuclei.

2. **Segmentation** (per FOV): Runs Cellpose in **2D+stitch mode** -- each
   z-slice is segmented independently as a 2D image (`do_3D=False`), then
   masks are stitched across z-slices using an IoU overlap threshold
   (`stitch_threshold`). If a cell in slice z=5 overlaps sufficiently with a
   cell in z=6, they receive the same cell ID. This is faster and more
   memory-efficient than true 3D segmentation while still producing
   consistent cell labels across the z-stack.

**Output:**

| File | Contents |
|------|----------|
| `segmentation/preprocessed/fov_*_preprocessed.tif` | Preprocessed 4-D volumes (nuclei + cytoplasm channels) |
| `segmentation/masks/fov_*_masks.tif` | Cell segmentation masks -- uint16 label images where each non-zero value = one cell |

**Config (3D mode, default):**

```yaml
segmentation:
  enabled: true
  aligned_images_dir: "/path/to/merlin_output/FiducialCorrelationWarp/images"  # required
  mode: "3d"                # segment all z-slices with 2D+stitch (default)
  nuclei_bit: 17            # 0-indexed bit for nuclei channel in dataorganization
  total_bits: 18            # total bits including nuclei + cell channels
  exclude_bits: []          # additional 0-indexed bits to exclude from cytoplasm
                            # (e.g. [0] to skip fiducial beads channel)
  median_kernel: 3          # median filter kernel size (odd integer)
  model_type: "cpsam"       # Cellpose model ("cyto2", "cpsam", etc.)
  diameter: null             # cell diameter in pixels (null = auto-detect)
  batch_size: 8             # GPU batch size
  stitch_threshold: 0.5     # IoU threshold for stitching 2D masks across z-slices
```

**Config (2D mode -- single z-slice):**

```yaml
segmentation:
  enabled: true
  aligned_images_dir: "/path/to/merlin_output/FiducialCorrelationWarp/images"
  mode: "2d"                # segment only one z-slice
  reference_z_slice: 5      # which z-slice to segment (1-indexed by default)
  z_indexing: 1              # 1 = 1-indexed (default), 0 = 0-indexed
  nuclei_bit: 17
  total_bits: 18
  exclude_bits: [0]         # e.g. skip fiducial beads (bit 0)
  model_type: "cpsam"
```

**Output format:**

- **3D mode (default):** The mask for each FOV is a 3-D array `(Z, Y, X)`
  -- one 2D label image per z-slice. The stitching step assigns consistent cell
  IDs across z-slices. For a FOV with 41 z-slices and 200 detected cells, the
  mask would be `(41, 2048, 2048)` uint16 with values 0-200.
- **2D mode:** Only the `reference_z_slice` is segmented. The output mask is a
  2-D array `(Y, X)`. Downstream cell assignment uses this single mask for
  barcodes from all z-slices. This is useful when only one z-slice has
  sufficient image quality for reliable segmentation (e.g. Nikon).

**Microscope compatibility:**

- **ONI / Nikon:** Tested and supported.
- **ANDOR:** Not yet validated. Current ANDOR image quality may not be
  sufficient for reliable cell segmentation. Compatibility will be assessed
  as ANDOR imaging protocols improve.

**Requires:** MERlin completed externally. Install with `pip install -e ".[segmentation]"`.
GPU is recommended but the stage falls back to CPU automatically if GPU is unavailable.

---

## Typical Workflows

### ONI or NIKON (standard)

```yaml
pipeline:
  stages: [index, stitch, focus_qc, inspect_positions, convert, merlin_config]
```

```bash
# 1. Run pre-MERlin stages
merfish-pipe run my_experiment.yaml

# 2. Run MERlin externally
source output/merlin_analysis/.merlinenv
bash output/merlin_analysis/run_merLIN.sh

# 3. Run post-MERlin analysis
merfish-pipe run my_experiment.yaml --stage correlation
```

### Andor

```yaml
pipeline:
  stages: [ims_convert, inspect_positions, merlin_config]
```

Andor uses `ims_convert` to convert IMS (HDF5) files to merged TIFFs and
extract stage positions. `inspect_positions` reads the per-round position
CSVs produced by `ims_convert` (no need to run `index` first). Stages
like `stitch`, `focus_qc`, `reregistration`, and `convert` are not used
in the Andor workflow.

### With reregistration

```yaml
pipeline:
  stages: [index, stitch, focus_qc, inspect_positions, reregistration, convert, merlin_config]

reregistration:
  enabled: true

filter_barcodes:
  enabled: true
```

```bash
# 1. Run pre-MERlin stages (including reregistration)
merfish-pipe run my_experiment.yaml

# 2. Run MERlin externally
source output/merlin_analysis/.merlinenv
bash output/merlin_analysis/run_merLIN.sh

# 3. Filter duplicate barcodes
merfish-pipe run my_experiment.yaml --stage filter_barcodes

# 4. Correlation (uses filtered barcodes automatically)
merfish-pipe run my_experiment.yaml --stage correlation
```

---

## Output Directory Layout

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
├── convert/                       (or ims_convert/)
│   └── merFISH_merged_*.tiff
├── parameters/                    (MERlin parameter files)
│   ├── dataorganization/
│   │   └── data_organization_{name}.csv
│   ├── microscope/
│   │   └── microscope_{name}.json
│   ├── positions/
│   │   └── positions_{name}.csv
│   ├── analysis/
│   │   └── analysis_{name}.json
│   └── codebooks/
│       └── {codebook}.csv
├── merlin_analysis/               (MERlin working directory)
│   ├── .merlinenv
│   ├── run_merLIN.sh
│   └── {experiment}/              (created by MERlin after it runs)
│       └── ExportBarcodes/barcodes.csv
├── merlin_config/                 (stage metadata only)
│   └── run_metadata.json
├── filter_barcodes/               (if enabled)
│   └── barcodes_filtered.csv
├── correlation/                   (if enabled)
│   ├── info_{name}.csv
│   └── correlation_plots.pdf
└── segmentation/                  (if enabled)
    ├── preprocessed/
    └── masks/
```

Every stage writes a `run_metadata.json` file in its output directory with
timestamp, status, runtime, and a list of output files.

---

## Troubleshooting

### Stage says "skipped" -- outputs already exist

The pipeline skips stages whose outputs are already present. To re-run:

```bash
merfish-pipe run my_experiment.yaml --stage <name> --force
```

### Config validation fails

Check the exact error message. Common causes:

- Unknown field name (typo in YAML key)
- Missing required field (`paths.raw_data_dir`, `paths.output_dir`, etc.)
- Invalid stage name in `pipeline.stages`

Use `merfish-pipe config show my_experiment.yaml` to see how your config was
resolved after merging.

### A stage fails

Check the stage's log output and `run_metadata.json`:

```bash
cat output_dir/<stage_name>/run_metadata.json
```

The `status` field will be `"failed"` with an `error` message. Fix the issue
and re-run with `--force`.

### "No images found" or "manifest is empty"

Your `paths.raw_data_dir` doesn't match the expected directory layout for your
microscope. Check:

- ONI/NIKON: expects channel subfolders with TIFF files matching the file pattern
- Andor: expects round-named folders with `.ims` files

### MERlin doesn't start

Make sure you sourced the environment file before running:

```bash
source output_dir/merlin_config/.merlinenv
```

Check that `run_merlin.sh` has the correct paths by reading it.

### trajectory_plot.html is missing

Install plotly: `pip install plotly` or `pip install -e ".[viz]"`. The pipeline
skips this plot if plotly is not installed.

### Segmentation fails with "no module named cellpose"

Install the segmentation extra: `pip install -e ".[segmentation]"`.

### cellpose install fails with "fastremap" build error

This happens on systems with old GCC (< 9.3). Use conda instead:
`conda install -c conda-forge cellpose`. Or install fastremap from a pre-built
wheel: `pip install --only-binary :all: fastremap` then `pip install cellpose`.

### Segmentation fails with "channel_axis and z_axis must be specified"

This pipeline requires cellpose v4.0+. Upgrade with:
`pip install --upgrade cellpose`.
