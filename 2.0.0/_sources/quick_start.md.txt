# Quick Start

This guide will help you get started with `pyremap` for remapping data between different spatial representations.

## Installation

To install `pyremap`, use `conda` (from the `conda-forge` channel):

```bash
conda install pyremap
```

## Basic Usage

Here is a simple example of how to use `pyremap` to remap data from a source grid to a destination grid.

### Step 1: Define Source and Destination Grids

```python
from pyremap import LatLonGridDescriptor, MpasCellMeshDescriptor

# Define a source grid (e.g., a regular latitude-longitude grid)
source_grid = LatLonGridDescriptor.read(file_name='source_grid.nc')

# Define a destination grid (e.g., an unstructured mesh)
destination_grid = MpasCellMeshDescriptor.read(file_name='destination_grid.nc')
```

### Step 2: Generate a Mapping File

```python
from pyremap import Remapper

# Create a remapper object
remapper = Remapper(map_filename='map_source_to_dest.nc', method='bilinear')
remapper.src_descriptor = source_grid
remapper.dst_descriptor = destination_grid

# Generate the mapping file
remapper.build_map()
```

### Step 3: Apply the Mapping

```python
# Remap data from the source grid to the destination grid
remapper.ncremap(
    in_filename='source_data.nc',
    out_filename='remapped_data.nc',
    variable_list=['temperature', 'salinity']
)
```

## Additional Resources

For more detailed examples and advanced usage, refer to the [documentation](index.md#pyremap).
