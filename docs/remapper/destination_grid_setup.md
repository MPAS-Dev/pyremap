# Destination Grid Setup

How to set up the destination grid for remapping.

## Overview
The destination grid defines the spatial representation of the output data after remapping. It can also be created from various grid types.

## Methods
- `dst_from_lon_lat`: Create a destination grid from latitude and longitude arrays.
- `dst_from_mpas`: Create a destination grid from an MPAS mesh file.
- `dst_from_points`: Create a destination grid from a collection of points.
- `dst_from_proj`: Create a destination grid from a projected grid.

## Example
```python
from pyremap import Remapper

remapper = Remapper()
remapper.dst_from_lon_lat(lat, lon)
```

## Setting `dst_descriptor` Directly
The `dst_descriptor` attribute can be set directly to define the destination
grid. This is useful when you already have a pre-defined grid descriptor
object or need functionality not provided by the `dst_from_*()` methods.

### Example
```python
from pyremap import Remapper, GridDescriptor

# Create or load a GridDescriptor object
dst_descriptor = GridDescriptor()

# Set the destination descriptor directly
remapper = Remapper()
remapper.dst_descriptor = dst_descriptor
```
