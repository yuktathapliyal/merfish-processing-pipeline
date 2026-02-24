# Reregistration & Data Conversion (Stages 5--8)

These stages transform your raw data into the format MERlin expects.
Reregistration (optional) corrects z-drift. The conversion stages merge raw
per-channel images into stacked multi-frame TIFFs. Finally, `merlin_config`
generates all the parameter files MERlin needs to run.

**After this group completes, you run MERlin externally.** The pipeline does
not run MERlin itself -- it prepares everything MERlin needs and gives you a
launch script.

---

## `reregistration`

**What it does:** When the microscope focus drifts during an experiment,
different FOVs end up with their "best focus" at different z-slices. This
means one FOV might have useful data at z=3..15, while another has it at
z=5..17. MERlin needs all FOVs to have the same z-range.

Reregistration fixes this by remapping z-slices so every FOV starts from the
same reference point. It reads the best-focus data from `focus_qc`, computes a
uniform target depth (the smallest usable z-range across all FOVs), and then
copies raw image files with new z-numbering into `remapped_data/`.

**When to use:** Only when `focus_qc` shows significant z-variation across
FOVs. Most common with ONI microscopes. Not needed for ANDOR.

**Output:**

| File | What it contains |
|------|-----------------|
| `reregistration/zmap_new_to_old.csv` | A table showing how each new z-index maps to the original z-index, per FOV. The `is_duplicate` column flags z-slices that had to be duplicated (copied from an adjacent slice) to fill gaps. |
| `remapped_data/{channel}/*.TIFF` | The remapped image files. Same format as the originals but with new z-numbering. These replace the raw images for all downstream stages. |

**What to check:** Open `run_metadata.json` and look at `target_z` -- this is
the uniform z-depth all FOVs were mapped to. In `zmap_new_to_old.csv`, check
how many slices are flagged as duplicates. A few duplicates are normal; if most
slices are duplicates, the original z-ranges were very inconsistent.

**Config:**

```yaml
reregistration:
  enabled: true         # must be explicitly enabled (default: false)
  total_z: null         # auto-detected from data (or override manually)
  target_z: null        # auto-computed from best-focus (or override manually)
```

**Requires:** `focus_qc` stage completed.

**Downstream effects:** When reregistration is enabled:
- `convert` automatically uses the remapped data instead of raw data
- `merlin_config` adjusts the data organization for the new z-count
- `filter_barcodes` should be run after MERlin to remove barcodes on
  duplicated z-slices

---

## `convert` (ONI and NIKON only)

**What it does:** ONI and NIKON microscopes save each channel and z-slice as a
separate TIFF file. MERlin expects a single stacked TIFF per (round, FOV) that
contains all channels and z-slices together. This stage merges those separate
files into the stacked format.

If reregistration was run, it automatically uses the remapped data from
`remapped_data/` instead of the original raw files.

**Output:**

| File | What it contains |
|------|-----------------|
| `merlin_data/merFISH_merged_{round}_{fov}.tiff` | One multi-frame stacked TIFF per (round, FOV). Each frame in the stack is one (channel, z-slice) combination. These go into `merlin_data/` which MERlin reads as its input directory. |

**What to check:** Verify the number of output TIFFs matches (rounds x FOVs).
You can open one in ImageJ and scroll through frames to confirm all channels
and z-slices are present.

**Config:** No special options needed.

**Requires:** `index` stage completed. Use `ims_convert` instead for ANDOR data.

---

## `ims_convert` (ANDOR only)

**What it does:** ANDOR microscopes save data as IMS (HDF5) files -- a
proprietary format that stores all channels and z-slices in a single
hierarchical file. This stage reads the IMS files, extracts the image data,
reorders channels if needed, and writes the same stacked TIFF format that
MERlin expects.

It also extracts stage position metadata from the IMS headers and writes
per-round position CSV files (used by `inspect_positions` and `merlin_config`).

**Output:**

| File | What it contains |
|------|-----------------|
| `merlin_data/merFISH_merged_{round}_{fov}.tiff` | Merged stacked TIFFs, same format as `convert` output. |
| `merlin_data/stagePos_Round#{round}.csv` | Stage position CSVs extracted from IMS metadata. One file per round, with x/y coordinates for each FOV. |

**Config:**

```yaml
raw_data:
  andor:
    channel_order: [0, 2, 1]    # reorder channels if your ANDOR setup needs it
```

**Requires:** Raw IMS files organised in round folders (e.g. `1st round/`,
`2nd round/`, `R1/`, `R2/`).

---

## `merlin_config`

**What it does:** Generates all the configuration files MERlin needs and
creates a shell script to launch MERlin. This is the last pipeline stage
before MERlin runs. It does NOT run MERlin itself -- just prepares everything.

MERlin needs several parameter files to know how to decode your experiment:
which bit corresponds to which image frame (data organization), microscope
optical parameters, FOV positions, analysis settings, and the codebook mapping
barcodes to genes.

**Output -- Parameter files** (saved to `parameters/`, organized by type):

| File | What it contains |
|------|-----------------|
| `parameters/dataorganization/data_organization_{name}.csv` | Maps each frame in the stacked TIFFs to a bit index. Tells MERlin which image frame corresponds to which imaging round, channel, and z-slice. |
| `parameters/microscope/microscope_{name}.json` | Microscope optical parameters: pixel size, image dimensions, flip/transpose flags. |
| `parameters/positions/positions_{name}.csv` | FOV stage positions (x, y) for MERlin's spatial coordinate system. |
| `parameters/analysis/analysis_{name}.json` | MERlin analysis task parameters: which algorithms to run and their settings. |
| `parameters/codebooks/{codebook}.csv` | The barcode codebook: maps each barcode ID to a gene name. Copied from your configured codebook template. |

**Output -- Launch files** (saved to `merlin_analysis/`):

| File | What it contains |
|------|-----------------|
| `merlin_analysis/.merlinenv` | Shell environment variables: sets `DATA_HOME`, `ANALYSIS_HOME`, and `PARAMETERS_HOME` so MERlin knows where to find everything. |
| `merlin_analysis/run_merLIN.sh` | A ready-to-run shell script that invokes MERlin with the correct arguments. |

**Config:**

```yaml
merlin:
  codebook_template: "/path/to/codebook.csv"   # required: barcode-to-gene mapping
  analysis_template: null                        # uses microscope default if null
  microscope_template: null                      # uses microscope default if null
  cores: 100                                     # number of cores for MERlin
```

**Requires:** `convert` or `ims_convert` stage completed.

> **Next step:** After `merlin_config` completes, run MERlin externally:
>
> ```bash
> source output_dir/merlin_analysis/.merlinenv
> bash output_dir/merlin_analysis/run_merLIN.sh
> ```
>
> MERlin will decode the barcodes and write its output to
> `merlin_analysis/{experiment}/ExportBarcodes/barcodes.csv`. Once MERlin
> finishes, continue with the [post-MERlin analysis stages](stages-post-merlin.md).

---

Next: [Post-MERlin Analysis Stages](stages-post-merlin.md)
