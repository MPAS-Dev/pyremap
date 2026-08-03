# Source Grid Setup
```{index} single: Remapper; Source Grid Setup
```

How to set up the source grid for remapping.

## Overview
The source grid defines the spatial representation of the data to be remapped. It can be created from various grid types, such as latitude-longitude grids, MPAS meshes, or projected grids.

## Methods
- `src_from_lon_lat`: Create a source grid from latitude and longitude arrays.
- `src_from_mpas`: Create a source grid from an MPAS mesh file.
- `src_from_proj`: Create a source grid from a projected grid.

## Example
```python
from pyremap import Remapper

remapper = Remapper()
remapper.src_from_lon_lat(lat, lon)
```

## Regional vs. Global Grids
By default, both `src_from_lon_lat` and `dst_from_lon_lat` determine whether a
1D lat/lon grid is regional or global automatically, based on whether the grid
is periodic in longitude (i.e. its longitudes wrap a full circle). Latitude
extent does not affect this: a zonally periodic but latitude-bounded grid (a
polar cap or zonal band) is treated as global so that longitude wraps across
the antimeridian.

The optional `regional` argument overrides this auto-detection. Pass
`regional=True` to force a grid to be treated as regional or `regional=False`
to force it to be treated as global:
```python
remapper.src_from_lon_lat('grid_file.nc', regional=False)
```

## Setting Source Descriptor Directly
The `src_descriptor` attribute can be set directly to define the source grid.
This is useful when you already have a pre-defined grid descriptor or need
additional capabilities not provided in the `src_from_*()` methods.

### Example
```python
from pyremap import Remapper, LatLonGridDescriptor

lat_lon_descriptor = LatLonGridDescriptor.read(file_name='grid_file.nc')
remapper = Remapper()
remapper.src_descriptor = lat_lon_descriptor
```

