# CLI Reference

The pipeline is controlled through the `merfish-pipe` command-line tool. It has
three main commands: `run`, `config`, and `status`.

---

## `merfish-pipe run`

Run pipeline stages.

```bash
merfish-pipe run <experiment.yaml> [OPTIONS]
```

### Flags

| Flag | What it does |
|------|-------------|
| `--stage <name>` | Run only this one stage (e.g. `--stage index`) |
| `--from-stage <name>` | Run from this stage onward, skipping earlier stages |
| `--dry-run` | Show what would run without actually executing anything |
| `--force` | Re-run stages even if their outputs already exist |
| `--profile slurm` | Submit stages as SLURM jobs instead of running locally |
| `--workers <N>` | Override the number of parallel worker threads |
| `--slurm` | Generate SLURM sbatch scripts with dependency chaining |
| `-v` | Enable verbose (DEBUG level) logging |

### Examples

```bash
# Run all stages listed in your config's pipeline.stages
merfish-pipe run my_experiment.yaml

# Run only the index stage
merfish-pipe run my_experiment.yaml --stage index

# Resume from stitch onward (useful after fixing an issue)
merfish-pipe run my_experiment.yaml --from-stage stitch

# Preview what would happen without running anything
merfish-pipe run my_experiment.yaml --dry-run

# Force re-run a stage (ignores existing output)
merfish-pipe run my_experiment.yaml --stage stitch --force

# Run on a SLURM cluster
merfish-pipe run my_experiment.yaml --profile slurm
```

---

## `merfish-pipe config`

Configuration utilities: generate templates, validate, and inspect configs.

### `config init`

Generate an experiment config template for your microscope:

```bash
# Minimal template (quick to fill in)
merfish-pipe config init --microscope oni -o my_experiment.yaml

# Detailed template (richly commented, explains every field and stage)
merfish-pipe config init --microscope oni --detailed -o my_experiment.yaml
```

| Flag | What it does |
|------|-------------|
| `--microscope` | Required. One of `oni`, `nikon`, `andor`. |
| `-o <path>` | Output file path. Prints to stdout if omitted. |
| `--detailed` | Include extensive comments explaining every field, stage, and common workflows. |

Without `--detailed`, you get a compact template with just the fields and
defaults -- good when you already know the config system. With `--detailed`,
you get the full reference template with inline documentation for every
option.

### `config validate`

Check your config for errors without running anything:

```bash
merfish-pipe config validate my_experiment.yaml
```

Prints experiment name, microscope type, paths, and stage list if valid.
Reports the exact error if something is wrong.

### `config show`

Display the fully resolved config (all three layers merged):

```bash
merfish-pipe config show my_experiment.yaml           # YAML format
merfish-pipe config show my_experiment.yaml --json    # JSON format
```

Useful for debugging -- see exactly what values the pipeline will use.

---

## `merfish-pipe status`

Check which stages have completed, failed, or are still pending:

```bash
merfish-pipe status my_experiment.yaml
```

Reads `run_metadata.json` from each stage's output directory and reports the
status. Shows timestamps and runtime for completed stages.

---

## Stage dependencies

Each stage needs certain inputs to exist before it can run. This table shows
what each stage requires and what it auto-detects from previous stages.

| Stage | Requires | Auto-detects from |
|-------|----------|-------------------|
| `index` | raw data directory | -- |
| `stitch` | `index` | -- |
| `focus_qc` | `index` | -- |
| `inspect_positions` | `index` OR `ims_convert` | -- |
| `reregistration` | `focus_qc` | -- |
| `convert` | `index` (+ `reregistration` if enabled) | `remapped_data_dir` |
| `ims_convert` | raw IMS files (ANDOR only) | -- |
| `merlin_config` | `convert` OR `ims_convert` | reregistration metadata |
| `filter_barcodes` | MERlin completed externally | `merlin_data_dir` |
| `correlation` | `filter_barcodes` OR MERlin + `bulk_file` | filtered barcodes if available |
| `optimize_correlation` | `correlation` | best threshold from `info_*.csv` |
| `joint_optimization` | `correlation` (current + all listed experiments) | best threshold from each experiment's `info_*.csv` |
| `segmentation` | MERlin completed + `aligned_images_dir` | -- |
| `cell_assignment` | `segmentation` + barcodes | masks dir, barcodes (filter_barcodes or MERlin) |
| `barcode_qc` | barcodes + codebook | barcodes: `cell_assignment` > `filter_barcodes` > MERlin |
| `anndata_export` | `cell_assignment` | -- |
| `spatial_visualization` | barcodes | barcodes: `cell_assignment` > `filter_barcodes` > MERlin |

> **Note:** MERlin runs externally (not a pipeline stage). After `merlin_config`
> generates the launch files, you run MERlin yourself, then resume the pipeline
> with `--from-stage filter_barcodes`.

### Running a single stage

```bash
# Run just barcode_qc (assumes all upstream stages / MERlin are done)
merfish-pipe run experiment.yaml --stage barcode_qc

# Run from correlation onwards (skips preprocessing)
merfish-pipe run experiment.yaml --from-stage correlation

# Dry-run a single stage to check that its inputs exist
merfish-pipe run experiment.yaml --stage segmentation --dry-run
```

### Auto-detection

Many post-MERlin stages auto-detect their input files so you don't need to
configure explicit paths. For example, `barcode_qc` checks for barcodes in
this order:

1. `cell_assignment/barcodes_assigned.csv` (preferred -- includes Cell_ID)
2. `filter_barcodes/barcodes_filtered.csv`
3. MERlin's `ExportBarcodes/barcodes.csv`

You can override auto-detection with explicit config (e.g.
`barcode_qc.barcodes_file: /path/to/barcodes.csv`).

---

Back to: [Documentation Index](README.md)
