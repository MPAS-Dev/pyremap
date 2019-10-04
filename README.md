# pyremap
[![Documentation Status](http://readthedocs.org/projects/pyremap/badge/?version=latest)](http://pyremap.readthedocs.io/en/latest/?badge=latest)

Python remapping tools for climate and earth system models.

## Documentation

[http://pyremap.readthedocs.io](http://pyremap.readthedocs.io)

## Installation

To use the latest version for developers, you will need to set up a conda
environment with the following packages:

 * python >= 3.6
 * numpy
 * scipy
 * netCDF4
 * xarray >= 0.10.0
 * dask
 * nco >= 4.8.1
 * pyproj

These can be installed via the conda command:
```
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n pyremap python=3.7 numpy scipy netCDF4 "xarray>=0.10.0" dask \
    "nco>=4.8.1" pyproj
conda activate pyremap
```

Then, get the code from:
[https://github.com/xylar/pyremap](https://github.com/xylar/pyremap)


## Examples

```
cd examples
```
The simplest way to make pyremap available is:
```
ln -s ../pyremap
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
