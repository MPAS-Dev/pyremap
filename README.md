# pyremap

[![Build Status](https://dev.azure.com/MPAS-Dev/pyremap%20testing/_apis/build/status/MPAS-Dev.pyremap?branchName=master)](https://dev.azure.com/MPAS-Dev/pyremap%20testing/_build/latest?definitionId=1&branchName=master)

Python remapping tools for climate and earth system models.

## Documentation

[http://mpas-dev.github.io/pyremap/stable/](http://mpas-dev.github.io/pyremap/stable/)

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
python -m pip install -e .
```

## Examples

```
cd examples
```
First, make mapping files for a lat-lon grid, and test it out by remapping
temperature from the example file:
```
wget https://web.lcrc.anl.gov/public/e3sm/inputdata/ocn/mpas-o/oQU240/ocean.QU.240km.151209.nc
./make_mpas_to_lat_lon_mapping.py
```
You should now see the mapping file:
```
map_oQU240_to_0.5x0.5degree.nc
```
as well as the input file (an initial condition for the MPAS-Ocean model) and
an example of temperature from the initial condition remapped to the new grid.
```
ocean.QU.240km.151209.nc
temp_0.5x0.5degree.nc
```

Second, let's try the same but to an Antarctic stereographic grid:
```
./make_mpas_to_Antarctic_stereo_mapping.py
```
Now, there's a new mapping file and example output file:
```
map_oQU240_to_6000.0x6000.0km_10.0km_Antarctic_stereo.nc
temp_6000.0x6000.0km_10.0km_Antarctic_stereo.nc
```

Finally, let's remap the temperature on the Antarctic grid to a lower
resolution Antarctic grid:
```
./remap_stereographic.py -i temp_6000.0x6000.0km_10.0km_Antarctic_stereo.nc \
    -o temp_6000.0x6000.0km_20.0km_Antarctic_stereo.nc -r 20
```
This created another mapping file and an output file:
```
map_6000x6000km_10km_Antarctic_stereo_to_6000x6000km_20.0km_Antarctic_stereo.nc
temp_6000.0x6000.0km_20.0km_Antarctic_stereo.nc
```
