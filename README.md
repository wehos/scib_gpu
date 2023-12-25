
# Benchmarking atlas-level data integration in single-cell genomics

This repository contains the code for a modified version of `scib` package.
For original instruction, refer to the [original repo](https://github.com/theislab/scib) and  [manuscript](https://doi.org/10.1038/s41592-021-01336-8).

# Environment Setup

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
