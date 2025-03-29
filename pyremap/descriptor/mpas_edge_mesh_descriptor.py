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


class MpasEdgeMeshDescriptor(MeshDescriptor):
    """
    A class for describing an MPAS edge mesh

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
            self.coords = {'latEdge': {'dims': 'nEdges',
                                       'data': ds.latEdge.values,
                                       'attrs': {'units': 'radians'}},
                           'lonEdge': {'dims': 'nEdges',
                                       'data': ds.lonEdge.values,
                                       'attrs': {'units': 'radians'}}}
            self.dims = ['nEdges']
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
        vertices_on_edge = ds_in.verticesOnEdge.values - 1
        cells_on_edge = ds_in.cellsOnEdge.values - 1
        dc_edge = ds_in.dcEdge.values
        dv_edge = ds_in.dvEdge.values
        sphere_radius = ds_in.attrs['sphere_radius']

        nedges = ds_in.sizes['nEdges']

        ds_out = xr.Dataset()

        valid_cells_on_edge = np.zeros(nedges)
        for icell in range(2):
            mask = cells_on_edge[:, icell] >= 0
            valid_cells_on_edge[mask] = valid_cells_on_edge[mask] + 1

        ds_out['grid_area'] = (
            ('grid_size',),
            0.5 * valid_cells_on_edge * dc_edge * dv_edge / (sphere_radius**2)
        )

        ds_out['grid_center_lat'] = (('grid_size',), lat_edge)
        ds_out['grid_center_lon'] = (('grid_size',), lon_edge)

        # grid corners:
        grid_corner_lon = np.zeros((nedges, 4))
        grid_corner_lat = np.zeros((nedges, 4))

        # start by repeating vertices, since these always exist
        vertices = vertices_on_edge[:, 0]
        grid_corner_lon[:, 0] = lon_vertex[vertices]
        grid_corner_lat[:, 0] = lat_vertex[vertices]
        grid_corner_lon[:, 1] = lon_vertex[vertices]
        grid_corner_lat[:, 1] = lat_vertex[vertices]
        vertices = vertices_on_edge[:, 1]
        grid_corner_lon[:, 2] = lon_vertex[vertices]
        grid_corner_lat[:, 2] = lat_vertex[vertices]
        grid_corner_lon[:, 3] = lon_vertex[vertices]
        grid_corner_lat[:, 3] = lat_vertex[vertices]

        # replace corners 0 and 2 with cells where they exist
        mask = cells_on_edge[:, 0] >= 0
        cells = cells_on_edge[mask, 0]
        grid_corner_lon[mask, 0] = lon_cell[cells]
        grid_corner_lat[mask, 0] = lat_cell[cells]
        mask = cells_on_edge[:, 1] >= 0
        cells = cells_on_edge[mask, 1]
        grid_corner_lon[mask, 2] = lon_cell[cells]
        grid_corner_lat[mask, 2] = lat_cell[cells]

        ds_out['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'), grid_corner_lat
        )
        ds_out['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'), grid_corner_lon
        )

        ds_out['grid_dims'] = xr.DataArray(
            [nedges],
            dims=('grid_rank',)
        ).astype('int32')

        ds_out['grid_imask'] = xr.DataArray(
            np.ones(nedges, dtype='int32'),
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
