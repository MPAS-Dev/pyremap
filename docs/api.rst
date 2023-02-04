#############
API reference
#############

This page provides an auto-generated summary of the pyremap API.

Mesh Descriptors
================

.. currentmodule:: pyremap.descriptor

.. autosummary::
   :toctree: generated/

   interp_extrap_corner
   interp_extrap_corners_2d

   MpasMeshDescriptor
   MpasMeshDescriptor.to_scrip
   MpasMeshDescriptor.to_esmf

   MpasEdgeMeshDescriptor
   MpasEdgeMeshDescriptor.to_scrip
   MpasEdgeMeshDescriptor.to_esmf

   get_lat_lon_descriptor
   LatLonGridDescriptor
   LatLonGridDescriptor.read
   LatLonGridDescriptor.create
   LatLonGridDescriptor.to_scrip
   LatLonGridDescriptor.to_esmf

   LatLon2DGridDescriptor
   LatLon2DGridDescriptor.read
   LatLon2DGridDescriptor.to_scrip

   ProjectionGridDescriptor
   ProjectionGridDescriptor.read
   ProjectionGridDescriptor.create
   ProjectionGridDescriptor.to_scrip
   ProjectionGridDescriptor.to_esmf

   PointCollectionDescriptor
   PointCollectionDescriptor.to_scrip
   PointCollectionDescriptor.to_esmf

Remapping
=========

.. currentmodule:: pyremap.remapper

.. autosummary::
   :toctree: generated/

   Remapper


Polar projections
=================

.. currentmodule:: pyremap.polar

.. autosummary::
   :toctree: generated/

   get_arctic_stereographic_projection
   get_antarctic_stereographic_projection
   get_polar_descriptor_from_file
   get_polar_descriptor
   to_polar
   from_polar

