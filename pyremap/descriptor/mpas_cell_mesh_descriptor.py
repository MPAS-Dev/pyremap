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
from pyremap.descriptor.utility import add_history, expand_scrip


class MpasCellMeshDescriptor(MeshDescriptor):
    """
    A class for describing an MPAS cell mesh

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
            self.coords = {
                'latCell': {
                    'dims': 'nCells',
                    'data': ds.latCell.values,
                    'attrs': {'units': 'radians'}
                },
                'lonCell': {
                    'dims': 'nCells',
                    'data': ds.lonCell.values,
                    'attrs': {'units': 'radians'}
                },
            }
            self.dims = ['nCells']
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
        vertices_on_cell = ds_in.verticesOnCell.values - 1
        nedges_on_cell = ds_in.nEdgesOnCell.values
        area_cell = ds_in.areaCell.values
        sphere_radius = ds_in.attrs['sphere_radius']

        ncells = ds_in.sizes['nCells']
        max_vertices = ds_in.sizes['maxEdges']

        ds_out = xr.Dataset()

        ds_out['grid_area'] = (
            ('grid_size',),
            area_cell / (sphere_radius**2)
        )

        ds_out['grid_center_lat'] = (('grid_size',), lat_cell)
        ds_out['grid_center_lon'] = (('grid_size',), lon_cell)

        # grid corners:
        grid_corner_lon = np.zeros((ncells, max_vertices))
        grid_corner_lat = np.zeros((ncells, max_vertices))
        for ivert in range(max_vertices):
            cell_indices = np.arange(ncells)
            # repeat the last vertex wherever iVertex > nEdgesOnCell
            local_vert_indices = np.minimum(nedges_on_cell - 1, ivert)
            vert_indices = vertices_on_cell[cell_indices, local_vert_indices]
            grid_corner_lat[cell_indices, ivert] = lat_vertex[vert_indices]
            grid_corner_lon[cell_indices, ivert] = lon_vertex[vert_indices]

        ds_out['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'), grid_corner_lat
        )
        ds_out['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'), grid_corner_lon
        )

        ds_out['grid_dims'] = xr.DataArray(
            [ncells],
            dims=('grid_rank',)
        ).astype('int32')

        ds_out['grid_imask'] = xr.DataArray(
            np.ones(ncells, dtype='int32'),
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
        self.write_netcdf(ds_out, scripFileName)
