# MpasVertexMeshDescriptor

A descriptor for MPAS vertex meshes.

## Overview
The `MpasVertexMeshDescriptor` class describes MPAS meshes based on vertex
geometry.

## Attributes
- `filename`: Path to the MPAS mesh file.
- `mesh_name`: Name of the MPAS mesh.

## Methods
- `to_scrip`: Converts the mesh to a SCRIP file.

## Example
```python
from pyremap import MpasVertexMeshDescriptor

descriptor = MpasVertexMeshDescriptor("mesh.nc", mesh_name="IcoswISC30E3r5")
descriptor.to_scrip("mesh.scrip.nc")
```