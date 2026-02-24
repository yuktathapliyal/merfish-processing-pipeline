# Configuration

## How the config system works

The pipeline uses a three-layer config system. Each layer provides settings,
and later layers override earlier ones:

```
Microscope defaults  -->  Your experiment YAML  -->  Execution profile  -->  Final config
   (built-in)              (you write this)           (built-in)
```

| Layer | What it controls | Where it lives | Do you edit it? |
|-------|-----------------|----------------|-----------------|
| **Microscope** | Pixel size, image dimensions, file patterns, orientation | `configs/microscopes/{oni,nikon,andor}.yaml` | No -- pre-built |
| **Experiment** | Your data paths, stage settings, codebook | Your YAML file (e.g. `my_experiment.yaml`) | **Yes** |
| **Execution** | Local vs SLURM, worker count | `configs/profiles/{local,slurm}.yaml` | No -- pre-built |

**In practice, you only need to create one file** -- your experiment config.
The microscope and execution defaults are already set up for you.

## Creating your experiment config

Generate an annotated template for your microscope type:

```bash
merfish-pipe config init --microscope oni -o my_experiment.yaml
```

Replace `oni` with `nikon` or `andor` to match your microscope. This creates a
YAML file pre-filled with sensible defaults and comments explaining each field.

## Key fields to fill in

Open the generated YAML file and edit these sections:

### Experiment identity

```yaml
experiment:
  name: "XP17596"              # a short unique ID for this experiment
  microscope: "oni"             # which microscope: oni, nikon, or andor
```

The `name` is used in output file names and MERlin parameters. Keep it short
and unique (e.g. your experiment number).

### Paths

```yaml
paths:
  raw_data_dir: "/data/raw/XP17596"      # where your raw images are
  output_dir: "/data/output/XP17596"     # where all pipeline results go
```

`raw_data_dir` should point to the top-level folder containing your raw images.
The pipeline will scan this directory to find images matching your microscope's
expected file patterns.

`output_dir` is where all pipeline outputs will be written. Each stage creates
its own subdirectory (e.g. `output_dir/index/`, `output_dir/stitch/`, etc.).

### Raw data settings

```yaml
raw_data:
  bead_channel_folder: "488nm, Raw"      # name of the fiducial channel subfolder
```

This tells the pipeline which channel subfolder contains the bead (fiducial)
images used for stitching and focus detection. The exact name depends on your
microscope and imaging setup.

### MERlin settings

```yaml
merlin:
  codebook_template: "/path/to/codebooks/C1E1_codebook.csv"
```

The codebook maps barcode IDs to gene names. This is only needed if you're
running `merlin_config` or the post-MERlin analysis stages.

### Pipeline stages

```yaml
pipeline:
  stages: [index, stitch, focus_qc, inspect_positions]
```

Lists which stages to run, in order. See the stage documentation for which
stages to include:
- [Pre-processing stages](stages-preprocessing.md)
- [Reregistration & conversion stages](stages-reregistration.md)
- [Post-MERlin stages](stages-post-merlin.md)

## Derived paths

These directories are auto-created under `output_dir`. You don't need to set
them unless you want to override the defaults:

| Directory | What it holds |
|-----------|--------------|
| `remapped_data/` | Reregistered images (only if reregistration is enabled) |
| `merlin_data/` | Merged TIFFs that MERlin will read as input |
| `merlin_analysis/` | MERlin launch files and MERlin's own output |
| `parameters/` | MERlin parameter files (data org, positions, codebook, etc.) |
| `logs/` | Pipeline log files |

## Validate your config

Before running anything, check that your config is valid:

```bash
merfish-pipe config validate my_experiment.yaml
```

If valid, this prints the experiment name, microscope type, paths, and stages.
If something is wrong, it tells you exactly what to fix.

## View the fully resolved config

To see what the final merged config looks like (all three layers combined):

```bash
merfish-pipe config show my_experiment.yaml
merfish-pipe config show my_experiment.yaml --json    # JSON format
```

This is useful for debugging -- you can see exactly what values the pipeline
will use, including defaults inherited from the microscope config.

---

Next: [Pre-processing Stages](stages-preprocessing.md)
