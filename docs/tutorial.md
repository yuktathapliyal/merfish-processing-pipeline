# merFISH Processing Pipeline -- Step-by-Step Tutorial

This tutorial walks you through installing, configuring, and running the
merFISH processing pipeline from scratch. It assumes no prior experience
with the tool.

> **Companion reference:** For a concise reference of all CLI flags, config
> fields, and output formats, see [user-manual.md](user-manual.md).

---

## Table of Contents

1. [Installation on HPC / Server](#installation-on-hpc--server)
2. [Phase 0: Verify Installation](#phase-0-verify-installation)
3. [Phase 1: Explore the Config System](#phase-1-explore-the-config-system)
4. [Phase 2: Create Your Experiment Config](#phase-2-create-your-experiment-config)
5. [Phase 3: Dry Run (Preview Without Executing)](#phase-3-dry-run-preview-without-executing)
6. [Phase 4: Run Stages One at a Time](#phase-4-run-stages-one-at-a-time)
7. [Phase 5: Run All Stages Together](#phase-5-run-all-stages-together)
8. [Phase 6: Post-MERlin Stages](#phase-6-post-merlin-stages)
9. [Quick Reference: What to Look For at Each Stage](#quick-reference-what-to-look-for-at-each-stage)
10. [Common Workflows](#common-workflows)

---

## Installation on HPC / Server

### Why conda?

Many HPC servers (e.g. Numbers at BC Cancer) run old system compilers like
GCC 4.8.5. Packages like `numpy`, `h5py`, and `opencv` need a modern compiler
to build from source. **Conda ships pre-built binaries** that don't need GCC
at all, so we use conda for the heavy compiled packages and pip only for our
pipeline package.

### Step 1: Create the conda environment

```bash
conda create -n merfish-pipe python=3.12 \
      numpy pandas h5py scikit-image opencv tifffile \
      matplotlib seaborn openpyxl pyyaml click tqdm \
      pytest pytest-cov \
      -c conda-forge -y
```

### Step 2: Activate and install pydantic via pip

```bash
conda activate merfish-pipe
pip install pydantic
```

### Step 3: Install the pipeline package

```bash
cd /path/to/merfish-processing-pipeline
pip install -e . --no-deps
```

The `--no-deps` flag tells pip to install **only our package** and trust that
conda already has everything. This avoids pip's dependency resolver entirely --
no more GCC failures.

### Optional extras

```bash
# Cell segmentation (requires GPU recommended)
pip install cellpose

# Interactive plotly trajectory plots
pip install plotly

# Development tools (if you plan to modify the code)
pip install ruff
```

### Local machine (alternative)

If you're on a machine with a modern compiler (macOS, Ubuntu, WSL2), you can
skip conda and use pip directly:

```bash
pip install -e ".[all,dev]"
```

---

## Phase 0: Verify Installation

After installation, run these three checks:

```bash
conda activate merfish-pipe

# 1. Version check
merfish-pipe --version

# 2. Help menu (should list: run, config, status)
merfish-pipe --help

# 3. Run the test suite
python -m pytest tests/ -v
```

All tests should pass. If you see import errors, the conda environment may be
missing a package -- check the error message and `conda install` what's needed.

---

## Phase 1: Explore the Config System

Before working with real data, get familiar with the config system.

### Generate templates for each microscope

```bash
merfish-pipe config init --microscope oni -o /tmp/test_oni.yaml
merfish-pipe config init --microscope nikon -o /tmp/test_nikon.yaml
merfish-pipe config init --microscope andor -o /tmp/test_andor.yaml
```

### Validate the templates

These will say VALID even though paths are placeholders -- validation checks
structure, not whether files exist on disk:

```bash
merfish-pipe config validate /tmp/test_oni.yaml
merfish-pipe config validate /tmp/test_nikon.yaml
merfish-pipe config validate /tmp/test_andor.yaml
```

### View the fully merged config

This shows all three layers (microscope + experiment + execution profile)
merged together:

```bash
merfish-pipe config show /tmp/test_oni.yaml
merfish-pipe config show /tmp/test_oni.yaml --json
```

### How the config system works

The pipeline uses three config layers that merge together (later wins):

| Layer | What it controls | Location |
|-------|-----------------|----------|
| **Microscope** | Pixel size, image dimensions, file patterns, flips | `configs/microscopes/{oni,nikon,andor}.yaml` |
| **Experiment** | Your data paths, stage settings, codebook | Your YAML file (you write this) |
| **Execution** | Local vs SLURM, worker count | `configs/profiles/{local,slurm}.yaml` |

You only write the experiment config. The microscope and execution configs
are pre-built.

---

## Phase 2: Create Your Experiment Config

Create a YAML file for your experiment. Here's a complete example:

```yaml
experiment:
  name: "XP14894"
  microscope: "nikon"         # "oni", "nikon", or "andor"

paths:
  raw_data_dir: "/data/raw/XP14894"        # where your raw images are
  output_dir: "/data/output/XP14894"       # where results go

raw_data:
  bead_channel_folder: "473nm, Raw"        # see table below

merlin:
  codebook_template: "/path/to/codebook.csv"  # needed for merlin_config stage

pipeline:
  stages: [index, stitch, focus_qc, inspect_positions]
  dry_run: false
  force: false
```

### Microscope-specific differences

| Setting | ONI | NIKON | ANDOR |
|---------|-----|-------|-------|
| `microscope:` | `"oni"` | `"nikon"` | `"andor"` |
| `bead_channel_folder:` | `"488nm, Raw"` | `"473nm, Raw"` | `"488nm, Raw"` |
| Position files | CSV (auto-detected) | XLSX (auto-detected) | Embedded in IMS |
| Convert stage | `convert` | `convert` | `ims_convert` |

### Validate your config

```bash
merfish-pipe config validate my_experiment.yaml
```

If valid you'll see:

```
VALID: experiment=XP14894 microscope=nikon
  output_dir: /data/output/XP14894
  raw_data_dir: /data/raw/XP14894
  execution: local (workers=8)
  stages: index, stitch, focus_qc, inspect_positions
```

---

## Phase 3: Dry Run (Preview Without Executing)

Before running anything for real, preview what each stage would do:

```bash
merfish-pipe run my_experiment.yaml --stage index --dry-run -v
merfish-pipe run my_experiment.yaml --stage stitch --dry-run -v
merfish-pipe run my_experiment.yaml --stage focus_qc --dry-run -v
merfish-pipe run my_experiment.yaml --stage inspect_positions --dry-run -v
merfish-pipe run my_experiment.yaml --stage merlin_config --dry-run -v
```

The `-v` flag enables verbose (DEBUG) logging so you can see exactly what the
stage is checking. Dry-run mode logs what would happen without writing output
files.

---

## Phase 4: Run Stages One at a Time

This is the core of the tutorial. Run each stage individually to understand
what it does and verify the outputs before moving on.

### Stage 1: `index` -- Scan Raw Data

**Why this stage:** Everything starts here. The pipeline needs to know what
files you have -- how many rounds, FOVs, z-slices, and channels. This stage
scans your raw data and writes a manifest.

```bash
merfish-pipe run my_experiment.yaml --stage index -v
```

**Check the outputs:**

```bash
ls output_dir/index/
head -5 output_dir/index/manifest.csv
head -5 output_dir/index/positions.standardized.csv
```

**What success looks like:**
- `manifest.csv` has columns: `round, fov, z_slice, channel, wavelength, abs_path, file_size`
- Each row is one raw image file
- `positions.standardized.csv` has columns: `round, tile_number, stage_pos_x, stage_pos_y`

**What failure looks like:**
- "No files found" -- your `raw_data_dir` doesn't match the expected layout
- Empty manifest -- wrong microscope type selected or wrong bead channel folder name

---

### Stage 2: `stitch` -- Build Tile Mosaics

**Why this stage:** Visualize how your FOVs tile together. The stitched
mosaics let you check tissue coverage, overlap quality, and whether the stage
positions are sensible.

```bash
merfish-pipe run my_experiment.yaml --stage stitch -v
```

**Check the outputs:**

```bash
ls output_dir/stitch/raw/
```

**What success looks like:**
- TIFF files like `IR01_mosaic.tiff`, `IR02_mosaic.tiff` (one per imaging round)
- Reasonable file sizes (a few MB each)
- Open in ImageJ/Fiji to visually inspect the mosaic

**What failure looks like:**
- Huge files (GB+) -- grid computation bug, tiles placed at wrong positions
- Missing files -- position data or bead images not found

**Depends on:** `index`

---

### Stage 3: `focus_qc` -- Best-Focus Detection

**Why this stage:** Each FOV has a z-stack. This finds which z-slice has the
sharpest focus in each FOV for each round. The results are used by
reregistration (if you enable it) and help you assess data quality.

```bash
merfish-pipe run my_experiment.yaml --stage focus_qc -v
```

**Check the outputs:**

```bash
head -5 output_dir/focus_qc/best_focus_slices.csv
# Open the heatmap to visually inspect z-variation:
# output_dir/focus_qc/heatmap.png
cat output_dir/focus_qc/summary.txt
```

**What success looks like:**
- `best_focus_slices.csv` has one row per FOV, one column per round (IR01, IR02, ...)
- Z values are in a valid range (e.g. 0-14 for a 15-slice stack)
- Heatmap shows gradual variation (not random noise)

**What failure looks like:**
- All zeros -- images couldn't be read or are blank
- "No TIFF files found" -- wrong bead channel folder

---

### Stage 4: `inspect_positions` -- Drift Analysis

**Why this stage:** The microscope stage moves between imaging rounds. This
measures how much each FOV drifted from its first-round position. If drift is
large, you may need reregistration.

```bash
merfish-pipe run my_experiment.yaml --stage inspect_positions -v
```

**Check the outputs:**

```bash
ls output_dir/inspect_positions/
head -5 output_dir/inspect_positions/drift_report.csv
cat output_dir/inspect_positions/drift_summary.txt
# Open drift_plot.png for visual inspection
# Open trajectory_plot.html in a browser (if plotly is installed)
```

**What success looks like:**
- `drift_report.csv` with per-FOV, per-round displacement values
- `drift_summary.txt` with aggregate statistics
- Drift plots showing consistent patterns

**What failure looks like:**
- Missing position or log file errors

**Depends on:** `index`

---

### Stage 5: `reregistration` -- Z-Depth Correction (Optional)

**Why this stage:** Different FOVs may have been imaged at different z-ranges
due to focus drift. Reregistration remaps z-slices so every FOV has the same
depth, making downstream analysis (MERlin) consistent.

**Enable it first** in your experiment YAML:

```yaml
reregistration:
  enabled: true
```

Then run it:

```bash
# Option A: Preview first (writes diagnostic CSVs without copying files)
merfish-pipe run my_experiment.yaml --stage reregistration --dry-run -v

# Option B: Run for real
merfish-pipe run my_experiment.yaml --stage reregistration -v

# Option C: Resume from reregistration onward (runs it + all later stages)
merfish-pipe run my_experiment.yaml --from-stage reregistration -v
```

**Check the outputs:**

```bash
head -10 output_dir/reregistration/zmap_new_to_old.csv
cat output_dir/reregistration/run_metadata.json   # check target_z value
ls output_dir/remapped_data/
```

**What success looks like:**
- `zmap_new_to_old.csv` with columns: FOV, IR, new_z, old_z, is_duplicate
- `remapped_data/` directory with channel subfolders containing remapped TIFFs
- File counts match expectations

**What failure looks like:**
- "focus CSV not found" -- `focus_qc` stage hasn't run yet
- Missing channel directories

**Depends on:** `focus_qc`

---

### Stage 6: `convert` -- Merge TIFFs for MERlin (ONI/NIKON)

**Why this stage:** Your raw data has separate TIFF files per channel, per
z-slice. MERlin expects a single stacked TIFF per (round, FOV) with all
channels and z-slices interleaved:

```
z0_wv0, z0_wv1, z0_wv2, z1_wv0, z1_wv1, z1_wv2, ...
```

```bash
merfish-pipe run my_experiment.yaml --stage convert -v
```

If reregistration was run, convert automatically reads from `remapped_data/`
instead of `raw_data_dir`.

**Check the outputs:**

```bash
ls output_dir/merlin_analysis/
# Should see: merFISH_merged_01_001.tiff, merFISH_merged_01_002.tiff, ...
```

**What success looks like:**
- One merged TIFF per (round, FOV) combination
- File sizes are reasonable (sum of all channels * z-slices)

**For ANDOR microscopes:** Use `ims_convert` instead of `convert`. It reads
IMS (HDF5) files and also extracts positions from the embedded metadata.

```bash
merfish-pipe run my_experiment.yaml --stage ims_convert -v
```

**Depends on:** `index` (+ `reregistration` if enabled)

---

### Stage 7: `merlin_config` -- Generate MERlin Parameter Files

**Why this stage:** MERlin needs several parameter files (data organization,
microscope parameters, analysis parameters, positions, codebook) and an
environment file. This stage generates all of them from your config and
templates.

```bash
merfish-pipe run my_experiment.yaml --stage merlin_config -v
```

**Check the outputs:**

```bash
# The run script that launches MERlin
cat output_dir/merlin_analysis/run_merLIN.sh

# The environment file
cat output_dir/merlin_analysis/.merlinenv

# Positions file (2-column, no header)
head -5 output_dir/parameters/positions/positions_XP14894.csv

# Data organization (expanded frame-to-bit mapping)
head -10 output_dir/parameters/dataorganization/data_organization_XP14894.csv
```

**What success looks like:**
- `run_merLIN.sh` uses the `-o` flag (not `-d`)
- `.merlinenv` sets DATA_HOME, ANALYSIS_HOME, PARAMETERS_HOME
- Positions CSV has exactly 2 columns (x, y) with no header
- If reregistration was used, DATA_HOME points to `remapped_data/` or `merlin_analysis/`

**What failure looks like:**
- Missing template errors (codebook, analysis, microscope JSONs)
- `-d` flag in run script (outdated MERlin version)
- Headers in positions CSV (MERlin expects no header)

**Depends on:** `convert` or `ims_convert`

---

### Running MERlin (External Step)

After `merlin_config` completes, run MERlin separately:

```bash
source output_dir/merlin_analysis/.merlinenv
bash output_dir/merlin_analysis/run_merLIN.sh
```

MERlin will write its outputs to `output_dir/merlin_analysis/{experiment_name}/`.
The key output for downstream stages is:

```
output_dir/merlin_analysis/{experiment_name}/ExportBarcodes/barcodes.csv
```

Wait for MERlin to complete before running post-MERlin stages.

---

## Phase 5: Run All Stages Together

Once you're comfortable with individual stages, run them all at once:

```bash
# Run all stages listed in pipeline.stages
merfish-pipe run my_experiment.yaml -v
```

### Force re-run

If outputs already exist, the pipeline skips those stages. To re-run:

```bash
# Force re-run all stages
merfish-pipe run my_experiment.yaml --force -v

# Force re-run just one stage
merfish-pipe run my_experiment.yaml --stage stitch --force -v
```

### Test the skip logic

```bash
# Run again without --force -- should skip everything
merfish-pipe run my_experiment.yaml -v
# Expected output: all stages show "skipping (outputs exist)"
```

### Check overall status

```bash
merfish-pipe status my_experiment.yaml
```

---

## Phase 6: Post-MERlin Stages

These stages run after MERlin has completed externally.

### `filter_barcodes` -- Remove Duplicate Z-Slice Barcodes

**Only needed if you used reregistration.** Reregistration pads some z-slices
by duplicating, which means MERlin may decode the same molecule twice. This
stage removes those duplicates.

Enable it in your config:

```yaml
filter_barcodes:
  enabled: true
  mode: "any"    # "any" = remove if ANY round duplicated, "all" = only if ALL rounds did
```

```bash
merfish-pipe run my_experiment.yaml --stage filter_barcodes -v
```

**Check:**

```bash
head -5 output_dir/filter_barcodes/barcodes_filtered.csv
cat output_dir/filter_barcodes/removal_summary.txt
```

---

### `correlation` -- Barcode vs Bulk RNA-seq Correlation

**Why:** Validates your MERFISH decoding by comparing barcode counts against
bulk RNA-seq expression data. Good Pearson/Spearman correlations (> 0.6)
indicate successful decoding.

Enable it in your config:

```yaml
correlation:
  enabled: true
  bulk_file: "/path/to/bulk_expression.csv"
```

```bash
merfish-pipe run my_experiment.yaml --stage correlation -v
```

**Check:**

```bash
cat output_dir/correlation/info_XP14894.csv
# Open output_dir/correlation/correlation_plots.pdf
```

---

### `segmentation` -- Cell Segmentation (Cellpose)

**Why:** Segments individual cells in the aligned microscopy images. Produces
per-FOV mask TIFFs where each pixel value is a cell ID.

Enable it in your config:

```yaml
segmentation:
  enabled: true
  aligned_images_dir: "/path/to/merlin_output/FiducialCorrelationWarp/images"
  nuclei_bit: 17
  total_bits: 18
```

```bash
merfish-pipe run my_experiment.yaml --stage segmentation -v
```

**Requires:** `pip install cellpose` (or `pip install -e ".[segmentation]"`)

---

## Quick Reference: What to Look For at Each Stage

| Stage | Success looks like | Failure looks like |
|-------|-------------------|-------------------|
| **index** | `manifest.csv` with round/fov/z/channel columns | "No files found" or empty manifest |
| **stitch** | `IR01_mosaic.tiff` etc., reasonable sizes (MBs) | Huge files (GB+) = grid bug |
| **focus_qc** | `best_focus_slices.csv` with z-values in valid range | All zeros or "no TIFF files" |
| **inspect_positions** | Drift report + plots | Missing position/log file |
| **reregistration** | `zmap_new_to_old.csv` + remapped TIFFs | "focus CSV not found" |
| **convert** | `merFISH_merged_*.tiff` in merlin_analysis/ | Missing manifest or source files |
| **ims_convert** | Same as convert, plus `stagePos_Round#*.csv` | IMS read errors |
| **merlin_config** | `run_merLIN.sh` with `-o` flag, 2-col positions CSV | Template not found, `-d` flag |
| **filter_barcodes** | `barcodes_filtered.csv` with fewer rows than input | Missing barcodes/zmap file |
| **correlation** | `info_*.csv` with Pearson r > 0.6 | Missing bulk expression file |
| **segmentation** | `masks/*.tiff` with cell labels (uint16 or int32) | "No module named cellpose" |

---

## Common Workflows

### ONI or NIKON (standard, no reregistration)

```yaml
pipeline:
  stages: [index, stitch, focus_qc, inspect_positions, convert, merlin_config]
```

```bash
# 1. Run pre-MERlin stages
merfish-pipe run my_experiment.yaml -v

# 2. Run MERlin externally
source output_dir/merlin_analysis/.merlinenv
bash output_dir/merlin_analysis/run_merLIN.sh

# 3. Run post-MERlin analysis
merfish-pipe run my_experiment.yaml --stage correlation -v
```

### ONI or NIKON (with reregistration)

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
merfish-pipe run my_experiment.yaml -v

# 2. Run MERlin externally
source output_dir/merlin_analysis/.merlinenv
bash output_dir/merlin_analysis/run_merLIN.sh

# 3. Filter duplicate barcodes from reregistration
merfish-pipe run my_experiment.yaml --stage filter_barcodes -v

# 4. Correlation (automatically uses filtered barcodes)
merfish-pipe run my_experiment.yaml --stage correlation -v
```

### ANDOR

```yaml
experiment:
  microscope: "andor"

pipeline:
  stages: [ims_convert, inspect_positions, merlin_config]

raw_data:
  bead_channel_folder: "488nm, Raw"
  andor:
    channel_order: [0, 2, 1]    # adjust for your Andor channel mapping
```

Same workflow as ONI/NIKON, but uses `ims_convert` instead of `convert`.
No `stitch` stage for Andor (IMS files handle tiling differently).

### The complete path to a MERlin-ready experiment

```bash
# Full pipeline (all stages in order)
merfish-pipe run my_experiment.yaml -v

# After MERlin completes, run post-analysis
merfish-pipe run my_experiment.yaml --stage filter_barcodes -v
merfish-pipe run my_experiment.yaml --stage correlation -v
merfish-pipe run my_experiment.yaml --stage segmentation -v
```

---

## Useful CLI Patterns

```bash
# Run a single stage
merfish-pipe run my_experiment.yaml --stage index

# Resume from a stage onward
merfish-pipe run my_experiment.yaml --from-stage convert

# Dry run (preview only)
merfish-pipe run my_experiment.yaml --dry-run

# Force re-run (ignore existing outputs)
merfish-pipe run my_experiment.yaml --stage stitch --force

# Verbose logging
merfish-pipe run my_experiment.yaml -v

# Check status
merfish-pipe status my_experiment.yaml

# Show fully resolved config
merfish-pipe config show my_experiment.yaml --json
```
