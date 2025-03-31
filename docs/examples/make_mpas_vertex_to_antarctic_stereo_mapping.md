# Make MPAS Vertex to Antarctic Stereo Mapping

:index:

This example demonstrates how to create a mapping file from MPAS mesh vertices to an Antarctic stereographic grid.

## Script Overview

The script uses `pyremap` to:
- Read an MPAS mesh file.
- Define a target Antarctic stereographic grid.
- Generate a mapping file for regridding based on mesh vertices.

## Example Code

```python
#!/usr/bin/env python

import numpy as np
import xarray as xr
from pyremap import Remapper, get_polar_descriptor

# Define the source MPAS mesh name and file
src_mesh_name = 'oQU240_vertex'
src_mesh_filename = 'ocean.QU.240km.151209.nc'

# Initialize the remapper with 4 tasks and the conservative regridding method
remapper = Remapper(ntasks=4, method='conserve')

# Set up the source mesh using MPAS vertex data
remapper.src_from_mpas(filename=src_mesh_filename, mesh_name=src_mesh_name, mesh_type='vertex')

# Define the target Antarctic stereographic grid descriptor
remapper.dst_descriptor = get_polar_descriptor(
    lx=6000.,  # Grid width in km
    ly=5000.,  # Grid height in km
    dx=10.,    # Grid cell size in x-direction (km)
    dy=10.,    # Grid cell size in y-direction (km)
    projection='antarctic'  # Antarctic stereographic projection
)

# Build the mapping file
remapper.build_map()
```

## Expected Output

The script generates a mapping file named:
```
map_oQU240_vertex_to_6000.0x5000.0km_10.0km_Antarctic_stereo_esmfaave.nc
```
that can be used for regridding MPAS vertex data to an Antarctic stereographic
grid.
