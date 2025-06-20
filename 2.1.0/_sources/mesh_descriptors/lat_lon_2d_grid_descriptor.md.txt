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

## Example
```python
from pyremap import LatLon2DGridDescriptor

descriptor = LatLon2DGridDescriptor.read("grid.nc")
descriptor.to_scrip("grid.scrip.nc")
```