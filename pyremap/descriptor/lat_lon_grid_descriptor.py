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
    interp_extrap_corner,
    round_res,
    unwrap_corners,
)


def get_lat_lon_descriptor(
    dlon, dlat, lon_min=-180.0, lon_max=180.0, lat_min=-90.0, lat_max=90.0
):
    """
    Get a descriptor of a lat-lon grid, used for remapping

    Parameters
    ----------
    dlon :  float
        Longitude resolution in degrees

    dlat :  float
        Latitude resolution in degrees

    lon_min :  float, optional
        Lower bound on longitude in degrees

    lon_max :  float, optional
        Upper bound on longitude in degrees

    lat_min :  float, optional
        Lower bound on latitude in degrees

    lat_max :  float, optional
        Upper bound on latitude in degrees

    Returns
    -------
    descriptor : pyremap.LatLonGridDescriptor object
        A descriptor of the lat/lon grid
    """
    nlat = int((lat_max - lat_min) / dlat) + 1
    nlon = int((lon_max - lon_min) / dlon) + 1
    lat = np.linspace(lat_min, lat_max, nlat)
    lon = np.linspace(lon_min, lon_max, nlon)

    descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

    return descriptor


class LatLonGridDescriptor(MeshDescriptor):
    """
    A class for describing a lat-lon grid

    Attributes
    ----------
    lat : Optional[numpy.ndarray]
        The latitude coordinate at grid-cell centers

    lon : Optional[numpy.ndarray]
        The longitude coordinate at grid-cell centers

    lat_corner : Optional[numpy.ndarray]
        The latitude coordinate at grid-cell corners

    lon_corner : Optional[numpy.ndarray]
        The longitude coordinate at grid-cell corners

    history : Optional[str]
        The history attribute written to SCRIP files
    """

    def __init__(self, mesh_name=None, regional=None):
        """
        Construct a mesh descriptor

        mesh_name : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid.  If ``None``, this will
            be determined automatically by checking the limits of the corner
            latitude and longitude to see if they cover the globe.
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
        regional=None,
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

        lat_var_name, lon_var_name : str, optional
            The name of the latitude and longitude variables in the grid file

        mesh_name : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid.  If ``None``, this will
            be determined automatically by checking the limits of the corner
            latitude and longitude to see if they cover the globe.
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
        descriptor.lon_corner = interp_extrap_corner(descriptor.lon)
        descriptor.lat_corner = interp_extrap_corner(descriptor.lat)

        descriptor._set_coords(
            lat_var_name,
            lon_var_name,
            ds[lat_var_name].dims[0],
            ds[lon_var_name].dims[0],
        )

        descriptor.history = add_history(ds=ds)
        return descriptor

    @classmethod
    def create(
        cls,
        lat_corner,
        lon_corner,
        units='degrees',
        mesh_name=None,
        regional=None,
    ):
        """
        Create the lat-lon grid with the given arrays and units.

        Parameters
        ----------
        lat_corner : numpy.ndarray
            One dimensional array defining the latitude coordinates of grid
            corners.

        lon_corner : numpy.ndarray
            One dimensional array defining the longitude coordinates of grid
            corners.

        units : {'degrees', 'radians'}, optional
            The units of `lat_corner` and `lon_corner`

        mesh_name : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid.  If ``None``, this will
            be determined automatically by checking the limits of the corner
            latitude and longitude to see if they cover the globe.
        """
        descriptor = cls(mesh_name=mesh_name, regional=regional)

        descriptor.lat_corner = lat_corner
        descriptor.lon_corner = lon_corner
        descriptor.lon = 0.5 * (lon_corner[0:-1] + lon_corner[1:])
        descriptor.lat = 0.5 * (lat_corner[0:-1] + lat_corner[1:])
        descriptor.units = units
        descriptor.history = add_history()
        descriptor._set_coords('lat', 'lon', 'lat', 'lon')
        return descriptor

    def to_scrip(self, scrip_filename, expand_dist=None, expand_factor=None):
        """
        Given a lat-lon grid file, create a SCRIP file based on the grid.

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
        assert self.lat_corner is not None, (
            'lat_corner must be set before calling to_scrip'
        )
        assert self.lon_corner is not None, (
            'lon_corner must be set before calling to_scrip'
        )
        assert self.units is not None, (
            'units must be set before calling to_scrip'
        )

        ds = xr.Dataset()

        (center_lon, center_lat) = np.meshgrid(self.lon, self.lat)
        (corner_lon, corner_lat) = np.meshgrid(
            self.lon_corner, self.lat_corner
        )

        ds['grid_center_lat'] = (('grid_size',), center_lat.flat)
        ds['grid_center_lon'] = (('grid_size',), center_lon.flat)
        ds['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'),
            unwrap_corners(corner_lat),
        )
        ds['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'),
            unwrap_corners(corner_lon),
        )

        ds['grid_dims'] = xr.DataArray(
            [len(self.lon), len(self.lat)], dims=('grid_rank',)
        ).astype('int32')

        ds['grid_imask'] = xr.DataArray(
            np.ones(ds.sizes['grid_size'], dtype='int32'), dims=('grid_size',)
        )

        ds.grid_center_lat.attrs['units'] = self.units
        ds.grid_center_lon.attrs['units'] = self.units
        ds.grid_corner_lat.attrs['units'] = self.units
        ds.grid_corner_lon.attrs['units'] = self.units
        ds.grid_imask.attrs['units'] = 'unitless'

        if expand_dist is not None or expand_factor is not None:
            expand_scrip(ds, expand_dist, expand_factor)

        ds.attrs['mesh_name'] = self.mesh_name
        ds.attrs['history'] = self.history
        self.write_netcdf(ds, scrip_filename)

    def _set_coords(
        self, lat_var_name, lon_var_name, lat_dim_name, lon_dim_name
    ):
        """
        Set up a coords dict with lat and lon
        """
        assert self.lat is not None, (
            'lat must be set before calling _set_coords'
        )
        assert self.lon is not None, (
            'lon must be set before calling _set_coords'
        )
        assert self.lat_corner is not None, (
            'lat_corner must be set before calling _set_coords'
        )
        assert self.lon_corner is not None, (
            'lon_corner must be set before calling _set_coords'
        )
        assert self.units is not None, (
            'units must be set before calling _set_coords'
        )

        self.lat_var_name = lat_var_name
        self.lon_var_name = lon_var_name
        self.coords = {
            lat_var_name: {
                'dims': lat_dim_name,
                'data': self.lat,
                'attrs': {'units': self.units},
            },
            lon_var_name: {
                'dims': lon_dim_name,
                'data': self.lon,
                'attrs': {'units': self.units},
            },
        }

        self.dims = [lat_dim_name, lon_dim_name]
        self.dim_sizes = [len(self.lat), len(self.lon)]

        # set the name of the grid
        dlat = self.lat[1] - self.lat[0]
        dlon = self.lon[1] - self.lon[0]
        lon_range = self.lon_corner[-1] - self.lon_corner[0]
        lat_range = self.lat_corner[-1] - self.lat_corner[0]
        if 'degree' in self.units:
            units = 'degree'
        elif 'rad' in self.units:
            units = 'radian'
        else:
            raise ValueError(f'Could not figure out units {self.units}')

        if self.regional is None:
            self.regional = False
            if units == 'degree':
                if np.abs(lon_range - 360.0) > 1e-10:
                    self.regional = True
                if np.abs(lat_range - 180.0) > 1e-10:
                    self.regional = True
            else:
                if np.abs(lon_range - 2.0 * np.pi) > 1e-10:
                    self.regional = True
                if np.abs(lat_range - np.pi) > 1e-10:
                    self.regional = True
        if self.mesh_name is None:
            self.mesh_name = (
                f'{round_res(abs(dlat))}x{round_res(abs(dlon))}{units}'
            )
