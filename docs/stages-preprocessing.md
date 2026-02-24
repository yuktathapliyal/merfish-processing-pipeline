# Pre-processing & QC (Stages 1--4)

These four stages scan your raw microscopy data and assess its quality. They
help you understand your data before making any modifications. Run these first
to check that images loaded correctly, tiles are positioned properly, focus is
consistent, and the microscope stage didn't drift too much between rounds.

**When to use:** Always. These are the starting point for every experiment.

---

## `index`

**What it does:** Scans your raw data directory and builds two reference files
that all downstream stages use. Think of this as the pipeline "learning" what
data you have -- how many rounds, FOVs, z-slices, and channels.

**Output:**

| File | What it contains |
|------|-----------------|
| `index/manifest.csv` | One row per raw image file. Columns: imaging round number, FOV number, z-slice number, channel name, and file path. This is the master inventory of your raw data. |
| `index/positions.standardized.csv` | Stage positions (x, y coordinates) for each FOV, in a unified format that works the same regardless of microscope type. Used by stitch and inspect_positions. |

**What to check:** Open `manifest.csv` and verify the number of unique rounds,
FOVs, and z-slices matches what you expect from your experiment. If rows are
missing, your `raw_data_dir` path or folder structure may be wrong.

**Config:** No special options. Just needs `paths.raw_data_dir` to point at
your raw images.

---

## `stitch`

**What it does:** Combines all FOV tiles from the bead (fiducial) channel into
a single large mosaic image. This gives you a bird's-eye view of your tissue
to visually check that all tiles are positioned correctly and the tissue is
where you expect it.

**Output:**

| File | What it contains |
|------|-----------------|
| `stitch/raw/*.TIFF` | One stitched mosaic per imaging round (or per z-slice, depending on `group_by`). Each file is a large TIFF showing all FOV tiles arranged in their correct spatial positions. |

If reregistration was run first, outputs go to `stitch/reregistered/` instead,
so you can compare stitches before and after z-correction.

**What to check:** Open the mosaics in ImageJ/FIJI. The tiles should fit
together without obvious gaps or overlaps. You should be able to see your
tissue sample across the mosaic. If tiles look scrambled, the position file
may be incorrect.

**Config:**

```yaml
stitch:
  group_by: "ir"        # "ir" = one mosaic per imaging round (default)
                        # "z"  = one mosaic per z-slice
  # ir_range: [1, 9]   # optional: only stitch these rounds
  # z_range: [1, 15]   # optional: only stitch these z-slices
```

**Requires:** `index` stage completed.

---

## `focus_qc`

**What it does:** For each FOV in each imaging round, scores every z-slice by
how "in focus" it is (using Laplacian-of-Gaussian variance on the bead channel)
and identifies the best-focus z-slice. This tells you at what depth the
microscope was best focused for each position.

If the best-focus z-slice varies a lot across FOVs, it means the microscope
focus drifted during the experiment -- this is when you'd want to enable
reregistration.

**Output:**

| File | What it contains |
|------|-----------------|
| `focus_qc/best_focus_slices.csv` | A table with columns: round, FOV, and best z-slice number. One row per (round, FOV) combination. |
| `focus_qc/heatmap.png` | A color-coded grid showing the best-focus z-slice for every FOV (rows) and round (columns). Uniform colors = good focus. A rainbow of colors = focus drift, consider reregistration. |
| `focus_qc/summary.txt` | Text summary: mean, std, and range of best-focus z-slices across all FOVs. |

**What to check:** Look at `heatmap.png`. If the colors are fairly uniform
(e.g. all green, meaning all FOVs focused at similar z-slices), your data is
well-focused. If you see a wide range of colors, the focus drifted and you
should consider running `reregistration`.

**Config:**

```yaml
focus_qc:
  sigma: 1.0    # Gaussian blur sigma (higher = smoother scoring)
  ksize: 3      # kernel size for the Laplacian filter
```

---

## `inspect_positions`

**What it does:** Compares the microscope stage positions across imaging rounds
to detect physical drift. If the microscope stage shifted between rounds, the
same FOV won't line up perfectly across rounds -- which can hurt barcode
decoding in MERlin.

**Output:**

| File | What it contains |
|------|-----------------|
| `inspect_positions/drift_report.csv` | Per-FOV, per-round position shift: delta_x (pixels), delta_y (pixels), and displacement (total movement in pixels). Round 1 is the reference. |
| `inspect_positions/drift_summary.txt` | Aggregate statistics: mean displacement, max displacement, worst FOV/round combination. |
| `inspect_positions/drift_plot.png` | Three-panel visualization: (1) displacement trend across rounds, (2) strip plot showing drift distribution per round, (3) heatmap of displacement by FOV and round. |
| `inspect_positions/trajectory_plot.html` | Interactive 3D plot showing the microscope stage path (requires `plotly`). Drag to rotate. Useful for spotting systematic drift patterns. |

**What to check:** Look at `drift_plot.png`. Small displacement values (< 1
pixel) are normal. Larger values mean the stage shifted significantly. The
heatmap panel shows whether drift is uniform across all FOVs or concentrated
in specific areas.

**Config:**

```yaml
inspect_positions:
  log_file: null                # auto-detected from raw data directory
  rounds_to_check: null         # check all rounds (or specify e.g. [1, 5, 9])
  trajectory_z_slices: null     # default: first 3 z-slices for trajectory plot
```

**Requires:** `index` stage completed, **or** `ims_convert` stage completed
(ANDOR workflow -- reads per-round `stagePos_Round#N.csv` files from
`merlin_data_dir`).

---

Next: [Reregistration & Conversion Stages](stages-reregistration.md)
