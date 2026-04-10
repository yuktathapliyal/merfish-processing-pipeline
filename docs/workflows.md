# Workflows

End-to-end examples for each microscope type. Pick the workflow that matches
your setup and follow the steps in order.

---

## ONI or NIKON (standard, no reregistration)

Use this workflow when `focus_qc` shows consistent focus across FOVs (no need
for z-correction).

**Config:**

```yaml
experiment:
  name: "XP17596"
  microscope: "oni"          # or "nikon"

paths:
  raw_data_dir: "/data/raw/XP17596"
  output_dir: "/data/output/XP17596"

raw_data:
  bead_channel_folder: "488nm, Raw"

merlin:
  codebook_template: "/path/to/codebooks/C1E1_codebook.csv"

pipeline:
  stages: [index, stitch, focus_qc, inspect_positions, convert, merlin_config]
```

**Steps:**

```bash
# 1. Run pre-MERlin stages
merfish-pipe run my_experiment.yaml

# 2. Check outputs: look at stitch mosaics, focus_qc heatmap, drift plots
#    If focus is inconsistent, consider the "with reregistration" workflow instead

# 3. Run MERlin externally
source output/merlin_analysis/.merlinenv
bash output/merlin_analysis/run_merLIN.sh

# 4. Run post-MERlin analysis (one stage at a time)
merfish-pipe run my_experiment.yaml --stage segmentation
merfish-pipe run my_experiment.yaml --stage cell_assignment
merfish-pipe run my_experiment.yaml --stage correlation
merfish-pipe run my_experiment.yaml --stage barcode_qc
merfish-pipe run my_experiment.yaml --stage anndata_export

# 5. Optional: find optimal gene subgroups and explore spatial data
merfish-pipe run my_experiment.yaml --stage optimize_correlation
merfish-pipe run my_experiment.yaml --stage joint_optimization   # requires joint_experiments in config
merfish-pipe run my_experiment.yaml --stage spatial_visualization
```

---

## ONI or NIKON (with reregistration)

Use this workflow when `focus_qc` shows significant z-variation across FOVs.
Reregistration corrects the z-drift before converting data for MERlin.

**Config:**

```yaml
experiment:
  name: "XP17596"
  microscope: "oni"

paths:
  raw_data_dir: "/data/raw/XP17596"
  output_dir: "/data/output/XP17596"

raw_data:
  bead_channel_folder: "488nm, Raw"

merlin:
  codebook_template: "/path/to/codebooks/C1E1_codebook.csv"

pipeline:
  stages: [index, stitch, focus_qc, inspect_positions, reregistration, convert, merlin_config]

reregistration:
  enabled: true

filter_barcodes:
  enabled: true

segmentation:
  enabled: true

cell_assignment:
  enabled: true
  crop_margin: 10             # filter cells near FOV edges
```

**Steps:**

```bash
# 1. Run pre-MERlin stages (including reregistration)
merfish-pipe run my_experiment.yaml

# 2. Check reregistration output: look at target_z in run_metadata.json
#    Optionally re-run stitch to compare before/after mosaics

# 3. Run MERlin externally
source output/merlin_analysis/.merlinenv
bash output/merlin_analysis/run_merLIN.sh

# 4. Filter duplicate barcodes from reregistration
merfish-pipe run my_experiment.yaml --stage filter_barcodes

# 5. Segment cells
merfish-pipe run my_experiment.yaml --stage segmentation

# 6. Assign barcodes to cells (uses filtered barcodes automatically)
merfish-pipe run my_experiment.yaml --stage cell_assignment

# 7. QC and analysis
merfish-pipe run my_experiment.yaml --stage correlation
merfish-pipe run my_experiment.yaml --stage barcode_qc
merfish-pipe run my_experiment.yaml --stage anndata_export

# 8. Optional: gene optimization and 3D visualization
merfish-pipe run my_experiment.yaml --stage optimize_correlation
merfish-pipe run my_experiment.yaml --stage joint_optimization   # requires joint_experiments in config
merfish-pipe run my_experiment.yaml --stage spatial_visualization
```

---

## ANDOR

ANDOR microscopes use a different file format (IMS/HDF5) and a simpler
workflow. No `index`, `stitch`, `focus_qc`, `reregistration`, or `convert`
stages are needed.

**Config:**

```yaml
experiment:
  name: "XP_ANDOR_001"
  microscope: "andor"

paths:
  raw_data_dir: "/data/raw/XP_ANDOR_001"
  output_dir: "/data/output/XP_ANDOR_001"

raw_data:
  bead_channel_folder: "488nm, Raw"
  andor:
    channel_order: [0, 2, 1]    # adjust for your ANDOR setup

merlin:
  codebook_template: "/path/to/codebooks/codebook.csv"

pipeline:
  stages: [ims_convert, inspect_positions, merlin_config]
```

**Steps:**

```bash
# 1. Convert IMS files and generate MERlin config
merfish-pipe run my_experiment.yaml

# 2. Check inspect_positions output for drift

# 3. Run MERlin externally
source output/merlin_analysis/.merlinenv
bash output/merlin_analysis/run_merLIN.sh

# 4. Post-MERlin analysis
merfish-pipe run my_experiment.yaml --stage correlation
merfish-pipe run my_experiment.yaml --stage barcode_qc

# 5. Optional: gene optimization and 3D visualization
merfish-pipe run my_experiment.yaml --stage optimize_correlation
merfish-pipe run my_experiment.yaml --stage joint_optimization   # requires joint_experiments in config
merfish-pipe run my_experiment.yaml --stage spatial_visualization
```

**Notes:**
- `ims_convert` handles both image conversion and position extraction
- `inspect_positions` reads the position CSVs produced by `ims_convert`
  (no need to run `index` first)
- Stages like `stitch`, `focus_qc`, `reregistration`, and `convert` are
  not used in the ANDOR workflow

---

Back to: [Documentation Index](README.md)
