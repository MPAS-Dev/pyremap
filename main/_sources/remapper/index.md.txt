# Remapper

:index:

Overview of the remapping process in pyremap.

Remapping is the process of interpolating data from one spatial representation
(source grid or mesh) to another (destination grid or mesh). This is a
critical step in climate and Earth system modeling, where data from different
models or observations often need to be compared or combined.

The remapping process in `pyremap` relies on **mesh descriptors** to define
the source and destination grid or mesh. A **mesh** can be unstructured,
such as those used in MPAS (Model for Prediction Across Scales), or
structured, such as regular latitude-longitude grids. Mesh descriptors provide
the metadata and functionality required to describe these grids or meshes,
enabling the creation of mapping (weight) files.

The remapping process typically involves the following steps:

1. **Source Grid Setup**:
   - Define the source grid or mesh using a mesh descriptor. This could be a
     latitude-longitude grid, an MPAS mesh, or a projected grid.
   - The source grid represents the spatial domain of the input data.

2. **Destination Grid Setup**:
   - Define the destination grid or mesh using a mesh descriptor. This could
     be a regular grid, a collection of points, or a polar stereographic
     projection.
   - The destination grid represents the spatial domain where the data will be
     interpolated.

3. **Optional Smoothing**:
   - If using the `conserve` mapping method, you can set the `expand_dist`
     and/or `expand_factor` attributes to define highly accurate and
     configurable smoothing to be applied as part of remapping.

4. **Building the Mapping**:
   - Generate a mapping (weight) file that defines how data from the source
     grid will be interpolated to the destination grid.
   - This step uses tools like ESMF (Earth System Modeling Framework) or MOAB
     to compute interpolation weights.

5. **Remapping Data**:
   - Apply the mapping file to remap data from the source grid to the
     destination grid.
   - This can be done using the NCO `ncremap` tool or through sparse matrix
     multiplication in `numpy`.

The `Remapper` class in `pyremap` provides a high-level interface for
performing these steps. It supports a variety of remapping methods, including
bilinear, conservative, and nearest-neighbor interpolation. The remapping
process is highly flexible, allowing users to customize the source and
destination grids, remapping methods, and tools.

```{toctree}
:maxdepth: 2

source_grid_setup
destination_grid_setup
smoothing
building_mapping
remapping_data
```
