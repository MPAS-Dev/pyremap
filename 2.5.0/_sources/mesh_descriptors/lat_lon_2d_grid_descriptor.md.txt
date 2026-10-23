# LatLon2DGridDescriptor
```{index} single: Mesh Descriptors; LatLon2DGridDescriptor
```

A descriptor for 2D latitude-longitude grids.

## Overview
The `LatLon2DGridDescriptor` class is used for grids where latitude and longitude are defined as 2D arrays.

## Attributes
- `lat_2d`: 2D array of latitude values.
- `lon_2d`: 2D array of longitude values.

## Methods
- `read`: Reads a 2D latitude-longitude grid from a file.
- `to_scrip`: Converts the grid to a SCRIP file.

## Grid-Cell Corners
As with {py:class}`LatLonGridDescriptor <pyremap.LatLonGridDescriptor>`,
`read()` uses the CF `bounds` of the latitude and longitude variables to find
grid-cell corners when they are available.  For 2D coordinates, the bounds
give the 4 vertices of each cell:
```
double lat(y, x) ;
    lat:units = "degrees_north" ;
    lat:bounds = "lat_bnds" ;
double lat_bnds(y, x, nv) ;
```
CF does not say which vertex comes first or which direction the 4 vertices are
traversed in, so pyremap works this out from the bounds themselves.  Both
latitude and longitude must have bounds, and neighboring cells must share
vertices, since the grid is described by 2D arrays of corners.  Otherwise,
corners are interpolated and extrapolated from cell centers and a warning is
raised.

## Example
```python
from pyremap import LatLon2DGridDescriptor

descriptor = LatLon2DGridDescriptor.read("grid.nc")
descriptor.to_scrip("grid.scrip.nc")
```