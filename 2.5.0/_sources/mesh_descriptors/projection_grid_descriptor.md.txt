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

## Grid-Cell Corners
`read()` uses the CF `bounds` of the `x` and `y` variables to find the corners
of each grid cell in projection space when they are available and describe
contiguous cells.  Otherwise, corners are interpolated between cell centers
and extrapolated at the ends of the grid.  Corners are transformed from
projection space to latitude and longitude by `to_scrip()`.

## Example
```python
from pyremap import ProjectionGridDescriptor

descriptor = ProjectionGridDescriptor.create(projection, x, y)
descriptor.to_scrip("grid.scrip.nc")
```