# ProjectionGridDescriptor
```{index} single: Mesh Descriptors; ProjectionGridDescriptor
```

A descriptor for grids based on map projections.

## Overview
The `ProjectionGridDescriptor` class describes grids defined by map projections.

## Attributes
- `projection`: The map projection used for the grid.

## Methods
- `read`: Reads a projection grid from a file.
- `create`: Creates a projection grid programmatically.
- `to_scrip`: Converts the grid to a SCRIP file.

## Example
```python
from pyremap import ProjectionGridDescriptor

descriptor = ProjectionGridDescriptor.create(projection, x, y)
descriptor.to_scrip("grid.scrip.nc")
```