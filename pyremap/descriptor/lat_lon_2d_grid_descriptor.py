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
from pyremap.descriptor.utility import (
    add_history,
    expand_scrip,
    interp_extrap_corners_2d,
    round_res,
    unwrap_corners,
)


class LatLon2DGridDescriptor(MeshDescriptor):
    """
    A class for describing a lat-lon grid that may not be a tensor grid
    (lat/lon are 2D arrays).  The grid is assumed to be regional, since this
    is difficult to determine just from the lat/lon values.  The calling code
    should pass ``regional=False`` to the constructor or ``read()`` method for
    global grids with 2D lat/lon.

    Attributes
    ----------
    lat : numpy.ndarray
        The latitude coordinate at grid-cell centers

    lon : numpy.ndarray
        The longitude coordinate at grid-cell centers

    lat_corner : numpy.ndarray
        The latitude coordinate at grid-cell corners

    lon_corner : numpy.ndarray
        The longitude coordinate at grid-cell corners

    history : str
        The history attribute written to SCRIP files
    """

    def __init__(self, mesh_name=None, regional=True):
        """
        Construct a mesh descriptor

        mesh_name : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid
        """
        super().__init__(mesh_name=mesh_name, regional=regional)
        self.lat: Optional[np.ndarray] = None
        self.lon: Optional[np.ndarray] = None
        self.units: Optional[str] = None
        self.lat_corner: Optional[np.ndarray] = None
        self.lon_corner: Optional[np.ndarray] = None
        self.history: Optional[str] = None

    @classmethod
    def read(
        cls,
        filename=None,
        ds=None,
        lat_var_name='lat',
        lon_var_name='lon',
        mesh_name=None,
        regional=True,
    ):
        """
        Read the lat-lon grid from a file with the given lat/lon var names.

        Parameters
        ----------
        filename : str, optional
            The path of the file containing the lat-lon grid (if ``ds`` is not
            supplied directly)

        ds : xarray.Dataset, optional
            The path of the file containing the lat-lon grid (if supplied,
            ``filename`` will be ignored)

        lat_var_name : str, optional
            The name of the latitude variable in the grid file

        lon_var_name : str, optional
            The name of the longitude variable in the grid file

        mesh_name : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid
        """
        if ds is None:
            ds = xr.open_dataset(filename)

        descriptor = cls(mesh_name=mesh_name, regional=regional)

        descriptor.mesh_name_from_attr(ds)
        # Get info from input file
        descriptor.lat = np.array(ds[lat_var_name].values, float)
        descriptor.lon = np.array(ds[lon_var_name].values, float)
        if 'degree' in ds[lat_var_name].units:
            descriptor.units = 'degrees'
        else:
            descriptor.units = 'radians'

        # interp/extrap corners
        descriptor.lon_corner = interp_extrap_corners_2d(descriptor.lon)
        descriptor.lat_corner = interp_extrap_corners_2d(descriptor.lat)

        descriptor._set_coords(
            lat_var_name,
            lon_var_name,
            ds[lat_var_name].dims[0],
            ds[lat_var_name].dims[1],
        )

        descriptor.history = add_history(ds=ds)
        return descriptor

    def to_scrip(self, scrip_filename, expand_dist=None, expand_factor=None):
        """
        Create a SCRIP file based on the grid.

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
        # Ensure required attributes are set
        assert self.lat is not None, 'lat must be set before calling to_scrip'
        assert self.lon is not None, 'lon must be set before calling to_scrip'
        assert self.lat_corner is not None, (
            'lat_corner must be set before calling to_scrip'
        )
        assert self.lon_corner is not None, (
            'lon_corner must be set before calling to_scrip'
        )
        assert self.units is not None, (
            'units must be set before calling to_scrip'
        )
        assert self.history is not None, (
            'history must be set before calling to_scrip'
        )

        ds = xr.Dataset()

        ds['grid_center_lat'] = (('grid_size',), self.lat.flat)
        ds['grid_center_lon'] = (('grid_size',), self.lon.flat)
        ds['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'),
            unwrap_corners(self.lat_corner),
        )
        ds['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'),
            unwrap_corners(self.lon_corner),
        )

        nlat, nlon = self.lat.shape

        ds['grid_dims'] = xr.DataArray(
            [nlon, nlat], dims=('grid_rank',)
        ).astype('int32')

        ds['grid_imask'] = xr.DataArray(
            np.ones(ds.sizes['grid_size'], dtype='int32'), dims=('grid_size',)
        )

        if expand_dist is not None or expand_factor is not None:
            expand_scrip(ds, expand_dist, expand_factor)

        ds.grid_center_lat.attrs['units'] = self.units
        ds.grid_center_lon.attrs['units'] = self.units
        ds.grid_corner_lat.attrs['units'] = self.units
        ds.grid_corner_lon.attrs['units'] = self.units
        ds.grid_imask.attrs['units'] = 'unitless'

        ds.attrs['mesh_name'] = self.mesh_name
        ds.attrs['history'] = self.history
        self.write_netcdf(ds, scrip_filename)

    def _set_coords(
        self, lat_var_name, lon_var_name, lat_dim_name, lon_dim_name
    ):
        """
        Set up a coords dict with lat and lon
        """
        # Ensure required attributes are set
        assert self.lat is not None, (
            'lat must be set before calling _set_coords'
        )
        assert self.lon is not None, (
            'lon must be set before calling _set_coords'
        )
        assert self.units is not None, (
            'units must be set before calling _set_coords'
        )

        self.lat_var_name = lat_var_name
        self.lon_var_name = lon_var_name
        self.coords = {
            lat_var_name: {
                'dims': (lat_dim_name, lon_dim_name),
                'data': self.lat,
                'attrs': {'units': self.units},
            },
            lon_var_name: {
                'dims': (lat_dim_name, lon_dim_name),
                'data': self.lon,
                'attrs': {'units': self.units},
            },
        }

        self.dims = [lat_dim_name, lon_dim_name]
        self.dim_sizes = self.lat.shape

        # set the name of the grid
        dlat = self.lat[1, 0] - self.lat[0, 0]
        dlon = self.lon[0, 1] - self.lon[0, 0]
        if 'degree' in self.units:
            units = 'degree'
        elif 'rad' in self.units:
            units = 'radian'
        else:
            raise ValueError(f'Could not figure out units {self.units}')
        if self.mesh_name is None:
            self.mesh_name = (
                f'{round_res(abs(dlat))}x{round_res(abs(dlon))}{units}'
            )
