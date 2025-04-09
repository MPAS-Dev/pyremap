# API reference
```{index} single: API reference
```

This page provides an auto-generated summary of the pyremap API.

## Mesh Descriptors

```{eval-rst}
.. currentmodule:: pyremap.descriptor

.. autosummary::
   :toctree: generated/


   interp_extrap_corner
   interp_extrap_corners_2d

   get_lat_lon_descriptor
   LatLonGridDescriptor
   LatLonGridDescriptor.read
   LatLonGridDescriptor.create
   LatLonGridDescriptor.to_scrip

   LatLon2DGridDescriptor
   LatLon2DGridDescriptor.read
   LatLon2DGridDescriptor.to_scrip

   MpasCellMeshDescriptor
   MpasCellMeshDescriptor.to_scrip

   MpasEdgeMeshDescriptor
   MpasEdgeMeshDescriptor.to_scrip

   MpasVertexMeshDescriptor
   MpasVertexMeshDescriptor.to_scrip

   PointCollectionDescriptor
   PointCollectionDescriptor.to_scrip

   ProjectionGridDescriptor
   ProjectionGridDescriptor.read
   ProjectionGridDescriptor.create
   ProjectionGridDescriptor.to_scrip
```

## Remapping

```{eval-rst}
.. currentmodule:: pyremap.remapper

.. autosummary::
   :toctree: generated/

   Remapper.src_from_lon_lat
   Remapper.src_from_mpas
   Remapper.src_from_proj


   Remapper.dst_from_lon_lat
   Remapper.dst_from_mpas
   Remapper.dst_from_points
   Remapper.dst_from_proj

   Remapper.build_map

   Remapper.ncremap
   Remapper.remap_numpy
```

## Polar projections

```{eval-rst}
.. currentmodule:: pyremap.polar

.. autosummary::
   :toctree: generated/

   get_arctic_stereographic_projection
   get_antarctic_stereographic_projection
   get_polar_descriptor_from_file
   get_polar_descriptor
   to_polar
   from_polar
```
