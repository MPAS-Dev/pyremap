# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2025 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2025 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2025 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/pyremap/main/LICENSE

import numpy as np
import xarray as xr

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import add_history, expand_scrip


class MpasEdgeMeshDescriptor(MeshDescriptor):
    """
    A class for describing an MPAS edge mesh

    Attributes
    ----------
    filename : str
        The path of the file containing the MPAS mesh

    history : str
        The history attribute written to SCRIP files
    """
    def __init__(self, filename, mesh_name=None):
        """
        Constructor stores the file name

        Parameters
        ----------
        filename : str
            The path of the file containing the MPAS mesh

        mesh_name : str, optional
            The name of the MPAS mesh (e.g. ``'oec60to30'`` or
            ``'orrs18to6'``).  If not provided, the data set in ``filename``
            must have a global attribute ``mesh_name`` that will be used
            instead.
        """
        super().__init__()

        with xr.open_dataset(filename) as ds:

            self.mesh_name = mesh_name
            self.mesh_name_from_attr(ds)
            if self.mesh_name is None:
                raise ValueError('No mesh_name provided or found in file.')

            self.filename = filename
            self.regional = True

            # build coords
            self.coords = {'lat_edge': {'dims': 'nEdges',
                                        'data': ds.latEdge.values,
                                        'attrs': {'units': 'radians'}},
                           'lon_edge': {'dims': 'nEdges',
                                        'data': ds.lonEdge.values,
                                        'attrs': {'units': 'radians'}}}
            self.dims = ['nEdges']
            self.dim_sizes = [ds.sizes[dim] for dim in self.dims]

            self.history = add_history(ds=ds)

    def to_scrip(self, scrip_filename, expand_dist=None, expand_factor=None):
        """
        Create a SCRIP file from the MPAS mesh.

        Parameters
        ----------
        scrip_filename : str
            The path to which the SCRIP file should be written

        expand_dist : float or numpy.ndarray, optional
            A distance in meters to expand each grid cell outward from the
            center.  If a ``numpy.ndarray``, one value per cell.

        expand_factor : float or numpy.ndarray, optional
            A factor by which to expand each grid cell outward from the center.
            If a ``numpy.ndarray``, one value per cell.
        """

        ds_in = xr.open_dataset(self.filename)
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

        n_edges = ds_in.sizes['nEdges']

        ds_out = xr.Dataset()

        valid_cells_on_edge = np.zeros(n_edges)
        for i_cell in range(2):
            mask = cells_on_edge[:, i_cell] >= 0
            valid_cells_on_edge[mask] = valid_cells_on_edge[mask] + 1

        ds_out['grid_area'] = (
            ('grid_size',),
            0.5 * valid_cells_on_edge * dc_edge * dv_edge / (sphere_radius**2)
        )

        ds_out['grid_center_lat'] = (('grid_size',), lat_edge)
        ds_out['grid_center_lon'] = (('grid_size',), lon_edge)

        # grid corners:
        grid_corner_lon = np.zeros((n_edges, 4))
        grid_corner_lat = np.zeros((n_edges, 4))

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
            [n_edges],
            dims=('grid_rank',)
        ).astype('int32')

        ds_out['grid_imask'] = xr.DataArray(
            np.ones(n_edges, dtype='int32'),
            dims=('grid_size',)
        )

        if expand_dist is not None or expand_factor is not None:
            expand_scrip(ds_out, expand_dist, expand_factor)

        ds_out.grid_center_lat.attrs['units'] = 'radians'
        ds_out.grid_center_lon.attrs['units'] = 'radians'
        ds_out.grid_corner_lat.attrs['units'] = 'radians'
        ds_out.grid_corner_lon.attrs['units'] = 'radians'
        ds_out.grid_imask.attrs['units'] = 'unitless'
        ds_out.grid_area.attrs['units'] = 'radians^2'

        ds_out.attrs['mesh_name'] = self.mesh_name
        ds_out.attrs['history'] = self.history
        self.write_netcdf(ds_out, scrip_filename)
