# Mesh Descriptors
```{index} single: Mesh Descriptors
```

Overview of mesh descriptors used in pyremap.

In order to perform remapping, we need a description of both the source and
destination mesh or grid. Here, we think of a mesh as being potentially
unstructured, whereas a grid is always structured. Mesh descriptors provide the
necessary metadata and functionality to describe these meshes or grids, enabling
remapping between different representations of spatial data.

The subclasses of `MeshDescriptor`, described in detail below, are used to
create grid files in
[SCRIP format](https://earthsystemmodeling.org/docs/release/ESMF_8_8_0/ESMF_refdoc/node3.html#SECTION03028100000000000000)
for use in generating mapping (weight) files with tools like
[ESMF](https://earthsystemmodeling.org/) or
[MOAB](https://sigma.mcs.anl.gov/moab-library/). These mapping files can then
be used to remap datasets using either the
[NCO `ncremap` tool](http://nco.sourceforge.net/nco.html#ncremap) or sparse
matrix multiplication in [`numpy`](https://numpy.org/).

```{toctree}
:maxdepth: 2

lat_lon_grid_descriptor
lat_lon_2d_grid_descriptor
mpas_cell_mesh_descriptor
mpas_edge_mesh_descriptor
mpas_vertex_mesh_descriptor
point_collection_descriptor
projection_grid_descriptor
```