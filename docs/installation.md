# Installation

## Core install

```bash
git clone <repo-url>
cd merfish-processing-pipeline
pip install -e .
```

This installs the pipeline with all core dependencies (numpy, pandas, tifffile,
opencv, matplotlib, etc.). No GPU or special hardware required for the core
stages.

## Optional extras

Some stages need additional packages. Install only what you need:

| Extra | What it adds | When you need it | Install command |
|-------|-------------|-----------------|----------------|
| `segmentation` | [Cellpose](https://github.com/MouseLand/cellpose) | Cell segmentation (stage 11) | `pip install -e ".[segmentation]"` |
| `export` | [AnnData](https://anndata.readthedocs.io/) | h5ad export for scanpy (stage 14) | `pip install -e ".[export]"` |
| `viz` | [Plotly](https://plotly.com/python/) | Interactive trajectory plots (stage 4) | `pip install -e ".[viz]"` |
| `all` | All of the above | Full pipeline | `pip install -e ".[all]"` |

## Cellpose GPU setup

Cellpose (used by the `segmentation` stage) runs **much** faster with a GPU.
If you have an NVIDIA GPU, install PyTorch with CUDA **before** the
segmentation extra:

```bash
# Check your CUDA version first
nvidia-smi

# Install PyTorch with CUDA (example for CUDA 12.1)
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121

# Then install the segmentation extra
pip install -e ".[segmentation]"

# Verify GPU is available
python -c "import torch; print(torch.cuda.is_available())"  # should print True
```

**CPU-only servers** (no GPU):

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install -e ".[segmentation]"
```

The pipeline will automatically fall back to CPU if no GPU is detected. It will
just be slower.

**If pip fails to build fastremap** (common on servers with old GCC < 9.3):

```bash
conda install -c conda-forge cellpose
```

## Server / HPC install (conda + no-deps)

On shared servers (e.g. Numbers, SLURM clusters), building C extensions from
source often fails because the system GCC is too old or you don't have root
access. The workaround is to let **conda** install the heavy compiled
dependencies, then install only the pipeline package with `--no-deps` so pip
doesn't try to rebuild anything.

```bash
# 1. Create a conda environment with all compiled dependencies
conda create -n merfish-pipe python=3.12 \
    numpy pandas h5py scikit-image opencv tifffile \
    matplotlib seaborn openpyxl pyyaml click tqdm \
    pytest pytest-cov \
    -c conda-forge -y

# 2. Activate the environment
conda activate merfish-pipe

# 3. Install pydantic via pip (not available on conda-forge with v2)
pip install pydantic

# 4. Install the pipeline itself -- no-deps tells pip to skip all
#    dependencies and trust that conda already has them
cd /path/to/merfish-processing-pipeline
pip install -e . --no-deps
```

The `--no-deps` flag is the key: it installs **only** our package and skips
pip's dependency resolver entirely. This avoids GCC build failures for
packages like fastremap, h5py, and scikit-image.

**Optional extras on the server:**

```bash
# Cellpose (for segmentation) -- install via conda to avoid build issues
conda install -c conda-forge cellpose -y

# AnnData (for h5ad export)
pip install anndata

# Plotly (for interactive trajectory plots)
pip install plotly
```

## Verify your install

```bash
merfish-pipe --help
```

You should see the list of available commands (`run`, `config`, `status`).

---

Next: [Configuration](configuration.md)
