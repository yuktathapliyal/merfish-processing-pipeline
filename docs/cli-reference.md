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

Generate an annotated experiment config template for your microscope:

```bash
merfish-pipe config init --microscope oni -o my_experiment.yaml
merfish-pipe config init --microscope nikon -o nikon_experiment.yaml
merfish-pipe config init --microscope andor -o andor_experiment.yaml
```

The template is pre-filled with sensible defaults and comments explaining each
field. Edit the required fields before running.

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

Back to: [Documentation Index](README.md)
