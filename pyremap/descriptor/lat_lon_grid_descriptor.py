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
from pyremap.descriptor.utility import (
    add_history,
    expand_scrip,
    interp_extrap_corner,
    round_res,
    unwrap_corners,
    write_netcdf,
)


def get_lat_lon_descriptor(dLon, dLat, lonMin=-180., lonMax=180., latMin=-90.,
                           latMax=90.):
    """
    Get a descriptor of a lat-lon grid, used for remapping

    Parameters
    ----------
    dLon :  float
        Longitude resolution in degrees

    dLat :  float
        Latitude resolution in degrees

    lonMin :  float, optional
        Lower bound on longitude in degrees

    lonMax :  float, optional
        Upper bound on longitude in degrees

    latMin :  float, optional
        Lower bound on latitude in degrees

    latMax :  float, optional
        Upper bound on latitude in degrees

    Returns
    -------
    descriptor : pyremap.LatLonGridDescriptor object
        A descriptor of the lat/lon grid
    """
    nLat = int((latMax - latMin) / dLat) + 1
    nLon = int((lonMax - lonMin) / dLon) + 1
    lat = np.linspace(latMin, latMax, nLat)
    lon = np.linspace(lonMin, lonMax, nLon)

    descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

    return descriptor


class LatLonGridDescriptor(MeshDescriptor):
    """
    A class for describing a lat-lon grid

    Attributes
    ----------
    lat : numpy.ndarray
        The latitude coordinate at grid-cell centers

    lon : numpy.ndarray
        The longitude coordinate at grid-cell centers

    latCorner : numpy.ndarray
        The latitude coordinate at grid-cell corners

    lonCorner : numpy.ndarray
        The longitude coordinate at grid-cell corners

    history : str
        The history attribute written to SCRIP files
    """
    def __init__(self, meshName=None, regional=None):
        """
        Construct a mesh descriptor

        meshName : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid.  If ``None``, this will
            be determined automatically by checking the limits of the corner
            latitude and longitude to see if they cover the globe.
        """
        super().__init__(meshName=meshName, regional=regional)
        self.lat = None
        self.lon = None
        self.units = None
        self.latCorner = None
        self.lonCorner = None
        self.history = None

    @classmethod
    def read(cls, fileName=None, ds=None, latVarName='lat',
             lonVarName='lon', meshName=None, regional=None):
        """
        Read the lat-lon grid from a file with the given lat/lon var names.

        Parameters
        ----------
        fileName : str, optional
            The path of the file containing the lat-lon grid (if ``ds`` is not
            supplied directly)

        ds : xarray.Dataset, optional
            The path of the file containing the lat-lon grid (if supplied,
            ``fileName`` will be ignored)

        latVarName, lonVarName : str, optional
            The name of the latitude and longitude variables in the grid file

        meshName : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid.  If ``None``, this will
            be determined automatically by checking the limits of the corner
            latitude and longitude to see if they cover the globe.
        """
        if ds is None:
            ds = xr.open_dataset(fileName)

        descriptor = cls(meshName=meshName, regional=regional)

        if descriptor.meshName is None and 'meshName' in ds.attrs:
            descriptor.meshName = ds.attrs['meshName']

        # Get info from input file
        descriptor.lat = np.array(ds[latVarName].values, float)
        descriptor.lon = np.array(ds[lonVarName].values, float)
        if 'degree' in ds[latVarName].units:
            descriptor.units = 'degrees'
        else:
            descriptor.units = 'radians'

        # interp/extrap corners
        descriptor.lonCorner = interp_extrap_corner(descriptor.lon)
        descriptor.latCorner = interp_extrap_corner(descriptor.lat)

        descriptor._set_coords(latVarName, lonVarName, ds[latVarName].dims[0],
                               ds[lonVarName].dims[0])

        descriptor.history = add_history(ds=ds)
        return descriptor

    @classmethod
    def create(cls, latCorner, lonCorner, units='degrees', meshName=None,
               regional=None):
        """
        Create the lat-lon grid with the given arrays and units.

        Parameters
        ----------
        latCorner : numpy.ndarray
            One dimensional array defining the latitude coordinates of grid
            corners.

        lonCorner : numpy.ndarray
            One dimensional array defining the longitude coordinates of grid
            corners.

        units : {'degrees', 'radians'}, optional
            The units of `latCorner` and `lonCorner`

        meshName : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid.  If ``None``, this will
            be determined automatically by checking the limits of the corner
            latitude and longitude to see if they cover the globe.
        """
        descriptor = cls(meshName=meshName, regional=regional)

        descriptor.latCorner = latCorner
        descriptor.lonCorner = lonCorner
        descriptor.lon = 0.5 * (lonCorner[0:-1] + lonCorner[1:])
        descriptor.lat = 0.5 * (latCorner[0:-1] + latCorner[1:])
        descriptor.units = units
        descriptor.history = add_history()
        descriptor._set_coords('lat', 'lon', 'lat', 'lon')
        return descriptor

    def to_scrip(self, scripFileName, expandDist=None, expandFactor=None):
        """
        Given a lat-lon grid file, create a SCRIP file based on the grid.

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
        ds = xr.Dataset()

        (center_lon, center_lat) = np.meshgrid(self.lon, self.lat)
        (corner_lon, corner_lat) = np.meshgrid(
            self.lonCorner, self.latCorner)

        ds['grid_center_lat'] = (('grid_size',), center_lat.flat)
        ds['grid_center_lon'] = (('grid_size',), center_lon.flat)
        ds['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'), unwrap_corners(corner_lat)
        )
        ds['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'), unwrap_corners(corner_lon)
        )

        ds['grid_dims'] = xr.DataArray(
            [len(self.lon), len(self.lat)],
            dims=('grid_rank',)
        ).astype('int32')

        ds['grid_imask'] = xr.DataArray(
            np.ones(ds.sizes['grid_size'], dtype='int32'),
            dims=('grid_size',)
        )

        if expandDist is not None or expandFactor is not None:
            expand_scrip(ds, expandDist, expandFactor)

        ds.grid_center_lat.attrs['units'] = self.units
        ds.grid_center_lon.attrs['units'] = self.units
        ds.grid_corner_lat.attrs['units'] = self.units
        ds.grid_corner_lon.attrs['units'] = self.units
        ds.grid_imask.attrs['units'] = 'unitless'

        ds.attrs['meshName'] = self.meshName
        ds.attrs['history'] = self.history
        write_netcdf(ds, scripFileName, format=self.format)

    def _set_coords(self, latVarName, lonVarName, latDimName,
                    lonDimName):
        """
        Set up a coords dict with lat and lon
        """
        self.latVarName = latVarName
        self.lonVarName = lonVarName
        self.coords = {latVarName: {'dims': latDimName,
                                    'data': self.lat,
                                    'attrs': {'units': self.units}},
                       lonVarName: {'dims': lonDimName,
                                    'data': self.lon,
                                    'attrs': {'units': self.units}}}

        self.dims = [latDimName, lonDimName]
        self.dimSize = [len(self.lat), len(self.lon)]

        # set the name of the grid
        dLat = self.lat[1] - self.lat[0]
        dLon = self.lon[1] - self.lon[0]
        lonRange = self.lonCorner[-1] - self.lonCorner[0]
        latRange = self.latCorner[-1] - self.latCorner[0]
        if 'degree' in self.units:
            units = 'degree'
        elif 'rad' in self.units:
            units = 'radian'
        else:
            raise ValueError('Could not figure out units {}'.format(
                self.units))

        if self.regional is None:
            self.regional = False
            if units == 'degree':
                if np.abs(lonRange - 360.) > 1e-10:
                    self.regional = True
                if np.abs(latRange - 180.) > 1e-10:
                    self.regional = True
            else:
                if np.abs(lonRange - 2. * np.pi) > 1e-10:
                    self.regional = True
                if np.abs(latRange - np.pi) > 1e-10:
                    self.regional = True
        if self.meshName is None:
            self.meshName = '{}x{}{}'.format(round_res(abs(dLat)),
                                             round_res(abs(dLon)), units)
