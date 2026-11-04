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
- `prepare_map`: Write the files the mapping tool needs and return the
  command that builds the mapping file, without running it.  This is for
  workflows that launch the tool themselves (see below).

## Attributes
- `map_filename`: The name of the file where the mapping weights will be saved.
- `method`: The remapping method to use (e.g., `bilinear`, `conserve`).
- `map_tool`: The tool to use for generating the mapping file (`esmf` or
  `moab`).
- `esmf_path`: The path to the ESMF installation, if `ESMF_RegridWeightGen` is
  not in your path.
- `moab_path`: The path to the MOAB installation, if `mptempest` is not in
  your path.
- `ntasks`: The number of MPI tasks to use to build the mapping file.  With
  `ntasks=1`, the mapping tool is run in serial (without `parallel_exec`) and,
  for `moab`, the meshes are not partitioned.
- `parallel_exec`: The parallel executable (e.g. `srun` or `mpirun`) to use
  to run the mapping tool.  This may include additional flags, e.g.
  `srun --mpi=pmi2`.
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

## Launching the mapping tool yourself
A workflow that runs MPI jobs with its own launcher and flags can call
`prepare_map()` instead of `build_map()`.  It writes the SCRIP files (and, for
`moab`, converts and partitions them) and returns the mapping tool's command
line.  The command has no parallel launcher in front of it, so the caller adds
one and runs it with `ntasks` MPI tasks.  `use_tmp` must be `False`, since the
tool reads the files after `prepare_map()` returns.

```python
import subprocess

from pyremap import Remapper

remapper = Remapper(ntasks=64, map_filename="map.nc", use_tmp=False)
remapper.map_tool = "moab"
remapper.src_from_lon_lat(lat, lon)
remapper.dst_from_lon_lat(lat_out, lon_out)
args = remapper.prepare_map()
subprocess.check_call(["srun", "-n", "64"] + args)
```

```{note}
For 1D lat/lon grids, whether a grid is treated as regional or global is
governed by longitude periodicity. Latitude extent does not force a grid to be
regional, so a zonally periodic but latitude-bounded grid is still global. Use
the `regional` argument of `src_from_lon_lat`/`dst_from_lon_lat` to override the
auto-detection.
```
