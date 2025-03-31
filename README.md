# pyremap

![pyremap logo](docs/_static/logo.png)

Python remapping tools for climate and earth system models.

## Documentation

[http://mpas-dev.github.io/pyremap/main/](http://mpas-dev.github.io/pyremap/main/)

## Installation

To use the latest version for developers, get the code from:
[https://github.com/MPAS-Dev/pyremap](https://github.com/MPAS-Dev/pyremap)

Then, you will need to set up a conda environment and install the package
in a way that points to the repo (so changes you make are available in the
conda environment):
```
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -y -n pyremap --file dev-spec.txt
conda activate pyremap
python -m pip install --no-deps --no-build-isolation -e .
```

## Examples

Detailed examples of how to use `pyremap` can be found in the [documentation](http://mpas-dev.github.io/pyremap/main/). These include:

1. Creating mapping files for remapping to a lat-lon grid.
2. Creating mapping files for remapping to an Antarctic stereographic grid.
3. Remapping data between different stereographic grids.

For step-by-step walkthroughs, see the [examples section](http://mpas-dev.github.io/pyremap/main/examples/index.html).
