# Remapping Data

How to remap data using the generated mapping file.

## Overview
Once the mapping file is created, it can be used to remap data from the source
grid to the destination grid.

## Methods
- `ncremap`: Use the NCO `ncremap` tool to remap data.
- `remap_numpy`: Use sparse matrix multiplication in `numpy` to remap data.

## Attributes
- `map_filename`: The mapping file to use for remapping data.
- `remap_tool`: The tool or method to use for remapping (e.g., `ncremap`, `remap_numpy`).

## Example
```python
from pyremap import Remapper

remapper = Remapper()
remapper.src_from_lon_lat(lat, lon)
remapper.dst_from_lon_lat(lat_out, lon_out)
remapper.build_map("map.nc", method="bilinear", map_tool="esmf")
remapper.remap_numpy("input.nc", "output.nc")
```
