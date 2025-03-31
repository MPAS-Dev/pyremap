# MpasEdgeMeshDescriptor

A descriptor for MPAS edge meshes.

## Overview
The `MpasEdgeMeshDescriptor` class describes MPAS meshes based on edge
geometry.

## Attributes
- `filename`: Path to the MPAS mesh file.
- `mesh_name`: Name of the MPAS mesh.

## Methods
- `to_scrip`: Converts the mesh to a SCRIP file.

## Example
```python
from pyremap import MpasEdgeMeshDescriptor

descriptor = MpasEdgeMeshDescriptor("mesh.nc", mesh_name="IcoswISC30E3r5")
descriptor.to_scrip("mesh.scrip.nc")
```