
# Benchmarking atlas-level data integration in single-cell genomics

This repository contains the code for a modified version of `scib` package.
For original instruction, refer to the [original repo](https://github.com/theislab/scib) and  [manuscript](https://doi.org/10.1038/s41592-021-01336-8).

# Benchmarking Workflow

## Step 1 Conda Installation
```
conda create -n leiden python=3.9
conda activate leiden
conda install cudatoolkit=11.7 -c pytorch -c nvidia
pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117
pip install --extra-index-url https://pypi.nvidia.com cudf-cu11==23.4.1 dask-cudf-cu11==23.4.1 cuml-cu11==23.4.1 cugraph-cu11==23.4.1 cucim==23.4.1
pip install einops ipdb pydance torchmetrics rapids-singlecell scvi-tools wandb hdf5plugin
conda install ipykernel
```

## Step 2 Scib Installation from Source

```
git clone git@github.com:wehos/scib_gpu.git
cd vllm
pip install -e .
```

## Step 3 Data Preparation
Prepare an `h5ad` file, with embeddings stored in `.obsm['X_emb']` or `.X`, batch labels in `batch_label`, cell type labels in `cell_type`.

## Step 4 Evaluation

Codes for standard evaluation.
```
INPUT_ANNDATA_PATH = 'xxx.h5ad'
OUTPUT_ANNDATA_PATH = 'yyy.csv'

import anndata as ad
import scib
import scanpy as sc

adata = ad.read_h5ad(INPUT_ANNDATA_PATH)
if 'X_emb' not in adata.obsm:
    adata.obsm['X_emb'] = adata.X      # If embeddings are stored in .X
if 'batch_label' not in adata.obs:
    adata.obs['batch_label'] = adata.obs['batch']

### To specify resolution choices for leiden, call `cluster_optimal_resolution` manually before scib.metrics
if False:    # Change to True to enable this code block
    sc.pp.neighbors(adata, use_rep='X_emb')
    scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", resolutions=[0.1, 0.25, 0.5, 0.75, 1.0], label_key="cell_type")
###

res = scib.metrics.metrics(adata, adata, "batch_label", "cell_type", embed="X_emb", cluster_key="cluster", organism='human',
                    ari_ = True, nmi_ = True, silhouette_=True, isolated_labels_asw_ =True, isolated_labels_=True, pcr_=True,
                    graph_conn_ =True, lisi_graph_=True)
res.to_csv(OUTPUT_ANNDATA_PATH)
```


# Instructions from original repo
### Please cite:

Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics.
Nat Methods 19, 41–50 (2022). [https://doi.org/10.1038/s41592-021-01336-8](https://doi.org/10.1038/s41592-021-01336-8)

## Package: scib

The `scib` python package is available on [PyPI](https://pypi.org/) and can be installed through

```commandline
pip install scib
```

Import `scib` in python:

```python
import scib
```

### Optional Dependencies

The package contains optional dependencies that need to be installed manually if needed.
These include R dependencies (`rpy2`, `anndata2ri`) which require an installation of R integration method packages.
All optional dependencies are listed under `setup.cfg` under `[options.extras_require]` and can be installed through pip.

e.g. for installing `rpy2` and `bbknn` dependencies:
```commandline
pip install 'scib[rpy2,bbknn]'
```

Optional dependencies outside of python need to be installed separately.
For instance, in order to run kBET, install it via the following command in R:

```R
install.packages('remotes')
remotes::install_github('theislab/kBET')
```

Install `pre-commit` to the repository for running it automatically every time you commit in git.

```shell
pre-commit install
```
