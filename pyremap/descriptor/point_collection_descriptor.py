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

from typing import Optional

import numpy as np
import xarray as xr

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import add_history


class PointCollectionDescriptor(MeshDescriptor):
    """
    A class for describing a collection of points

    Attributes
    ----------
    lat : Optional[numpy.ndarray]
        The latitude of each point

    lon : Optional[numpy.ndarray]
        The longitude of each point

    units : Optional[{'degrees', 'radians'}]
        The units of ``lats`` and ``lons``

    history : Optional[str]
        The history attribute written to SCRIP files
    """

    def __init__(
        self,
        lats,
        lons,
        collection_name,
        units='degrees',
        out_dimension='n_points',
    ):
        """
        Constructor stores

        Parameters
        ----------
        lats : numpy.ndarray
            The latitude of each point

        lons : numpy.ndarray
            The longitude of each point

        collection_name : str
            A unique name for the collection of transects, used in the names
            of files containing data mapped to these points.

        units : {'degrees', 'radians'}, optional
            The units of ``lats`` and ``lons``

        out_dimension : str, optional
            The name of the dimension corresponding to the points (i.e. the
            "horizontal" dimension of the point collection)
        """
        super().__init__(mesh_name=collection_name, regional=True)

        self.lat: Optional[np.ndarray] = lats
        self.lon: Optional[np.ndarray] = lons
        self.units: Optional[str] = units

        # build coords
        self.coords = {
            'lat': {
                'dims': out_dimension,
                'data': self.lat,
                'attrs': {'units': units},
            },
            'lon': {
                'dims': out_dimension,
                'data': self.lon,
                'attrs': {'units': units},
            },
        }
        self.dims = [out_dimension]
        self.dim_sizes = [len(self.lat)]
        self.history: Optional[str] = add_history()

    def to_scrip(self, scrip_filename, expand_dist=None, expand_factor=None):
        """
        Write a SCRIP file for the point collection

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

        assert self.lat is not None, 'lat must be set before calling to_scrip'
        assert self.lon is not None, 'lon must be set before calling to_scrip'
        assert self.units is not None, (
            'units must be set before calling to_scrip'
        )
        assert self.history is not None, (
            'history must be set before calling to_scrip'
        )

        ds = xr.Dataset()

        ds['grid_center_lat'] = (('grid_size',), self.lat)
        ds['grid_center_lon'] = (('grid_size',), self.lon)

        npoints = len(self.lat)
        # grid corners:
        grid_corner_lon = np.zeros((npoints, 4))
        grid_corner_lat = np.zeros((npoints, 4))
        # just repeat the center lat and lon
        for ivert in range(4):
            grid_corner_lat[:, ivert] = self.lat
            grid_corner_lon[:, ivert] = self.lon

        ds['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'),
            grid_corner_lat,
        )
        ds['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'),
            grid_corner_lon,
        )

        ds['grid_dims'] = xr.DataArray([npoints], dims=('grid_rank',)).astype(
            'int32'
        )

        ds['grid_imask'] = xr.DataArray(
            np.ones(npoints, dtype='int32'), dims=('grid_size',)
        )

        ds['grid_area'] = xr.DataArray(np.zeros(npoints), dims=('grid_size',))

        ds['grid_dims'] = (
            ('grid_rank',),
            [
                npoints,
            ],
        )
        ds['grid_imask'] = (('grid_size',), np.ones(npoints, dtype=int))

        ds.grid_center_lat.attrs['units'] = self.units
        ds.grid_center_lon.attrs['units'] = self.units
        ds.grid_corner_lat.attrs['units'] = self.units
        ds.grid_corner_lon.attrs['units'] = self.units
        ds.grid_imask.attrs['units'] = 'unitless'

        ds.attrs['mesh_name'] = self.mesh_name
        ds.attrs['history'] = self.history
        self.write_netcdf(ds, scrip_filename)
