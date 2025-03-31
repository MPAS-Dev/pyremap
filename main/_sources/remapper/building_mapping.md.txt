# Building the Mapping
```{index} single: Remapper; Building the Mapping
```

How to build the mapping file for remapping.

## Overview
The mapping file contains the weights needed to remap data from the source
grid to the destination grid. It can be generated using tools like ESMF or
MOAB.

## Methods
- `build_map`: Build the mapping file using the specified remapping method and
tool.

## Attributes
- `map_filename`: The name of the file where the mapping weights will be saved.
- `method`: The remapping method to use (e.g., `bilinear`, `conserve`).
- `map_tool`: The tool to use for generating the mapping file (`esmf` or
  `moab`).
- `esmf_path`: The path to the ESMF installation, if `ESMF_RegridWeightGen` is
  not in your path.
- `moab_path`: The path to the MOAB installation, if `mptempest` is not in
  your path.
- `parallel_exec`: The parallel executable (e.g. `srun` or `mpirun`) to use
  to run the mapping tool.
- `use_tmp`: Whether to use a temporary directory for intermediate files.

## Example
```python
from pyremap import Remapper

remapper = Remapper(map_filename="map.nc", method="bilinear")
remapper.map_tool = "esmf"
remapper.src_from_lon_lat(lat, lon)
remapper.dst_from_lon_lat(lat_out, lon_out)
remapper.build_map()
```
