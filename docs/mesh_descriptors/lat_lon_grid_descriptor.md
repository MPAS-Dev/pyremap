# LatLonGridDescriptor

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

## Example
```python
from pyremap import LatLonGridDescriptor

descriptor = LatLonGridDescriptor.create(lat, lon)
descriptor.to_scrip("grid.scrip.nc")
```