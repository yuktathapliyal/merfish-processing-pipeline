# Troubleshooting

Common errors and how to fix them.

---

## Stage says "skipped" -- outputs already exist

The pipeline automatically skips stages whose outputs already exist. To force
a re-run:

```bash
merfish-pipe run my_experiment.yaml --stage <name> --force
```

---

## Config validation fails

Check the exact error message. Common causes:

- **Unknown field name** -- typo in a YAML key (e.g. `stagees` instead of `stages`)
- **Missing required field** -- forgot to set `paths.raw_data_dir`, `paths.output_dir`, etc.
- **Invalid stage name** -- misspelled stage in `pipeline.stages`
- **Extra field** -- added a field the config doesn't recognize (the pipeline
  uses strict validation)

Use `merfish-pipe config show my_experiment.yaml` to see how your config was
resolved after merging all three layers.

---

## A stage fails

Check the stage's `run_metadata.json`:

```bash
cat output_dir/<stage_name>/run_metadata.json
```

The `status` field will be `"failed"` and the `error` field will contain the
error message. Fix the underlying issue and re-run with `--force`:

```bash
merfish-pipe run my_experiment.yaml --stage <name> --force
```

---

## "No images found" or "manifest is empty"

Your `paths.raw_data_dir` doesn't match the expected directory layout.

- **ONI/NIKON:** Expects channel subfolders (e.g. `488nm, Raw/`, `647nm, Raw/`)
  containing TIFF files matching the file pattern.
- **ANDOR:** Expects round-named folders (e.g. `1st round/`, `2nd round/` or
  `R1/`, `R2/`) containing `.ims` files.

Check that `raw_data_dir` points to the correct top-level folder.

---

## MERlin doesn't start

Make sure you sourced the environment file **before** running the launch script:

```bash
source output_dir/merlin_analysis/.merlinenv
bash output_dir/merlin_analysis/run_merLIN.sh
```

If it still fails, read `run_merLIN.sh` to check that the paths inside are
correct.

---

## trajectory_plot.html is missing

The interactive trajectory plot requires plotly. Install it:

```bash
pip install plotly
# or
pip install -e ".[viz]"
```

The pipeline skips this plot silently if plotly is not installed. All other
`inspect_positions` outputs still work.

---

## Segmentation fails with "no module named cellpose"

Cellpose is an optional dependency. Install the segmentation extra:

```bash
pip install -e ".[segmentation]"
```

---

## cellpose install fails with "fastremap" build error

This happens on systems with old GCC (< 9.3). Solutions:

1. Use conda: `conda install -c conda-forge cellpose`
2. Install fastremap from a pre-built wheel first:
   `pip install --only-binary :all: fastremap` then `pip install cellpose`

---

## Segmentation fails with "channel_axis and z_axis must be specified"

This pipeline requires cellpose v4.0 or later. Upgrade:

```bash
pip install --upgrade cellpose
```

---

## anndata_export writes CSV but no h5ad file

The h5ad export requires the `anndata` package. Install the export extra:

```bash
pip install -e ".[export]"
```

The CSV fallback (`cell_gene_matrix.csv`) always works without anndata.

---

Back to: [Documentation Index](README.md)
