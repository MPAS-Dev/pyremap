# Make MPAS to Lat-Lon Mapping
```{index} single: Examples; MPAS to Lat-Lon
```

This example demonstrates how to create a mapping file from an MPAS cell mesh
to a regular latitude-longitude grid.

## Script Overview

The script uses `pyremap` to:
- Read an MPAS mesh file.
- Define a target latitude-longitude grid.
- Generate a mapping file for regridding.

## Example Code

```python
from pyremap import Remapper, get_lat_lon_descriptor

# Input MPAS mesh file and name
src_mesh_filename = 'ocean.QU.240km.151209.nc'
src_mesh_name = 'oQU240'

# Configure the remapper
remapper = Remapper(ntasks=1, method='bilinear')
remapper.src_from_mpas(filename=src_mesh_filename, mesh_name=src_mesh_name)

# Define the target lat-lon grid
remapper.dst_descriptor = get_lat_lon_descriptor(dlon=0.5, dlat=0.5)

# Build the mapping file
remapper.build_map()

# Perform remapping (example with Python remapping)
src_filename = f'temp_{src_mesh_name}.nc'
out_filename = f'temp_{remapper.dst_descriptor.mesh_name}.nc'
remapper.ncremap(src_filename, out_filename)
```

## Expected Output

The script generates a mapping file `map_oQU240_to_0.5x0.5degree_esmfbilin.nc`
and a regridded output file `temp_0.5x0.5degree.nc` on the latitude-longitude
grid.
