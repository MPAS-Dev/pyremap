# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/pyremap/main/LICENSE

import numpy as np
import xarray as xr

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import add_history, expand_scrip, write_netcdf


class MpasVertexMeshDescriptor(MeshDescriptor):
    """
    A class for describing an MPAS vertex mesh

    Attributes
    ----------
    fileName : str
        The path of the file containing the MPAS mesh

    history : str
        The history attribute written to SCRIP files
    """
    def __init__(self, fileName, meshName=None):
        """
        Constructor stores the file name

        Parameters
        ----------
        fileName : str
            The path of the file containing the MPAS mesh

        meshName : str, optional
            The name of the MPAS mesh (e.g. ``'oEC60to30'`` or
            ``'oRRS18to6'``).  If not provided, the data set in ``fileName``
            must have a global attribute ``meshName`` that will be used
            instead.
        """
        super().__init__()

        with xr.open_dataset(fileName) as ds:

            if meshName is None:
                if 'meshName' not in ds.attrs:
                    raise ValueError('No meshName provided or found in file.')
                self.meshName = ds.attrs['meshName']
            else:
                self.meshName = meshName

            self.fileName = fileName
            self.regional = True

            # build coords
            self.coords = {'latVertex': {'dims': 'nVertices',
                                         'data': ds.latVertex.values,
                                         'attrs': {'units': 'radians'}},
                           'lonVertex': {'dims': 'nVertices',
                                         'data': ds.lonVertex.values,
                                         'attrs': {'units': 'radians'}}}
            self.dims = ['nVertices']
            self.dimSize = [ds.sizes[dim] for dim in self.dims]

            self.history = add_history(ds=ds)

    def to_scrip(self, scripFileName, expandDist=None, expandFactor=None):
        """
        Create a SCRIP file from the MPAS mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written

        expandDist : float or numpy.ndarray, optional
            A distance in meters to expand each grid cell outward from the
            center.  If a ``numpy.ndarray``, one value per cell.

        expandFactor : float or numpy.ndarray, optional
            A factor by which to expand each grid cell outward from the center.
            If a ``numpy.ndarray``, one value per cell.
        """

        ds_in = xr.open_dataset(self.fileName)
        lat_cell = ds_in.latCell.values
        lon_cell = ds_in.lonCell.values
        lat_vertex = ds_in.latVertex.values
        lon_vertex = ds_in.lonVertex.values
        lat_edge = ds_in.latEdge.values
        lon_edge = ds_in.lonEdge.values
        cells_on_vertex = ds_in.cellsOnVertex.values - 1
        edges_on_vertex = ds_in.edgesOnVertex.values - 1
        kite_areas_on_vertex = ds_in.kiteAreasOnVertex.values
        nvertices = ds_in.sizes['nVertices']
        vertex_degree = ds_in.sizes['vertexDegree']
        sphere_radius = ds_in.attrs['sphere_radius']

        ds_out = xr.Dataset()

        if vertex_degree != 3:
            raise ValueError(f'MpasVertexMeshDescriptor does not support '
                             f'vertexDegree {vertex_degree}')

        valid_cells_on_vertex = np.zeros(nvertices, dtype=int)
        vertex_area = np.zeros(nvertices)
        for icell in range(vertex_degree):
            mask = cells_on_vertex[:, icell] >= 0
            valid_cells_on_vertex[mask] = valid_cells_on_vertex[mask] + 1
            vertex_area[mask] = (
                vertex_area[mask] + kite_areas_on_vertex[mask, icell]
            )

        ds_out['grid_area'] = (
            ('grid_size',),
            vertex_area / (sphere_radius**2)
        )

        ds_out['grid_center_lat'] = (('grid_size',), lat_vertex)
        ds_out['grid_center_lon'] = (('grid_size',), lon_vertex)

        # grid corners:
        grid_corner_lon = np.zeros((nvertices, 6))
        grid_corner_lat = np.zeros((nvertices, 6))

        # start by filling the corners with the vertex lat and lon as the
        # fallback
        for icorner in range(6):
            grid_corner_lon[:, icorner] = lon_vertex
            grid_corner_lat[:, icorner] = lat_vertex

        # if edges on vertex exist, replace even indices with their locations
        for iedge in range(3):
            mask = edges_on_vertex[:, iedge] >= 0
            edges = edges_on_vertex[mask, iedge]
            grid_corner_lon[mask, 2 * iedge] = lon_edge[edges]
            grid_corner_lat[mask, 2 * iedge] = lat_edge[edges]

        # similarly, if cells on vertex exist, replace odd indices with their
        # locations
        for icell in range(3):
            mask = cells_on_vertex[:, icell] >= 0
            cells = cells_on_vertex[mask, icell]
            grid_corner_lon[mask, 2 * icell + 1] = lon_cell[cells]
            grid_corner_lat[mask, 2 * icell + 1] = lat_cell[cells]

        ds_out['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'), grid_corner_lat
        )
        ds_out['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'), grid_corner_lon
        )

        ds_out['grid_dims'] = xr.DataArray(
            [nvertices],
            dims=('grid_rank',)
        ).astype('int32')

        ds_out['grid_imask'] = xr.DataArray(
            np.ones(nvertices, dtype='int32'),
            dims=('grid_size',)
        )

        if expandDist is not None or expandFactor is not None:
            expand_scrip(ds_out, expandDist, expandFactor)

        ds_out.grid_center_lat.attrs['units'] = 'radians'
        ds_out.grid_center_lon.attrs['units'] = 'radians'
        ds_out.grid_corner_lat.attrs['units'] = 'radians'
        ds_out.grid_corner_lon.attrs['units'] = 'radians'
        ds_out.grid_imask.attrs['units'] = 'unitless'
        ds_out.grid_area.attrs['units'] = 'radians^2'

        ds_out.attrs['meshName'] = self.meshName
        ds_out.attrs['history'] = self.history
        write_netcdf(ds_out, scripFileName, format=self.format)
