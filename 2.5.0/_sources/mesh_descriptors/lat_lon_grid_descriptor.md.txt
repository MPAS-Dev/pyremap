# LatLonGridDescriptor
```{index} single: Mesh Descriptors; LatLonGridDescriptor
```

A descriptor for latitude-longitude grids.

## Overview
The `LatLonGridDescriptor` class is used to describe a regular latitude-longitude grid.

## Attributes
- `lat`: Latitude values.
- `lon`: Longitude values.

## Methods
- `read`: Reads a latitude-longitude grid from a file.
- `create`: Creates a latitude-longitude grid programmatically.
- `to_scrip`: Converts the grid to a SCRIP file.

## Grid-Cell Corners
Remapping needs the corners of each grid cell, not just the cell centers.
When `read()` is used, corners come from the CF `bounds` attribute of the
latitude and longitude variables if it is present:
```
double lat(lat) ;
    lat:units = "degrees_north" ;
    lat:bounds = "lat_bnds" ;
double lat_bnds(lat, nbnd) ;
```
The bounds must describe contiguous cells (the upper edge of each cell is the
lower edge of the next), since a 1D lat/lon grid is described by 1D arrays of
corners.  If the `bounds` attribute is missing, points to a variable that is
not in the file, has the wrong shape, or describes cells with gaps or overlaps
between them, corners are instead interpolated between cell centers and
extrapolated at the ends of the grid, and a warning is raised in all but the
first of these cases.

## Example
```python
from pyremap import LatLonGridDescriptor

descriptor = LatLonGridDescriptor.create(lat, lon)
descriptor.to_scrip("grid.scrip.nc")
```