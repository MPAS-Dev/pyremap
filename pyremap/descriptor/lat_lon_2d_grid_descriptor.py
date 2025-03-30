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

    latCorner : numpy.ndarray
        The latitude coordinate at grid-cell corners

    lonCorner : numpy.ndarray
        The longitude coordinate at grid-cell corners

    history : str
        The history attribute written to SCRIP files
    """
    def __init__(self, meshName=None, regional=True):
        """
        Construct a mesh descriptor

        meshName : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid
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
             lonVarName='lon', meshName=None, regional=True):
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

        latVarName : str, optional
            The name of the latitude variable in the grid file

        lonVarName : str, optional
            The name of the longitude variable in the grid file

        meshName : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names

        regional : bool or None, optional
            Whether this is a regional or global grid
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
        descriptor.lonCorner = interp_extrap_corners_2d(descriptor.lon)
        descriptor.latCorner = interp_extrap_corners_2d(descriptor.lat)

        descriptor._set_coords(latVarName, lonVarName, ds[latVarName].dims[0],
                               ds[latVarName].dims[1])

        descriptor.history = add_history(ds=ds)
        return descriptor

    def to_scrip(self, scripFileName, expandDist=None, expandFactor=None):
        """
        Create a SCRIP file based on the grid.

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

        ds['grid_center_lat'] = (('grid_size',), self.lat.flat)
        ds['grid_center_lon'] = (('grid_size',), self.lon.flat)
        ds['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'), unwrap_corners(self.latCorner)
        )
        ds['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'), unwrap_corners(self.lonCorner)
        )

        nlat, nlon = self.lat.shape

        ds['grid_dims'] = xr.DataArray(
            [nlon, nlat],
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
        self.write_netcdf(ds, scripFileName)

    def _set_coords(self, latVarName, lonVarName, latDimName,
                    lonDimName):
        """
        Set up a coords dict with lat and lon
        """
        self.latVarName = latVarName
        self.lonVarName = lonVarName
        self.coords = {latVarName: {'dims': (latDimName, lonDimName),
                                    'data': self.lat,
                                    'attrs': {'units': self.units}},
                       lonVarName: {'dims': (latDimName, lonDimName),
                                    'data': self.lon,
                                    'attrs': {'units': self.units}}}

        self.dims = [latDimName, lonDimName]
        self.dimSize = self.lat.shape

        # set the name of the grid
        dLat = self.lat[1, 0] - self.lat[0, 0]
        dLon = self.lon[0, 1] - self.lon[0, 0]
        if 'degree' in self.units:
            units = 'degree'
        elif 'rad' in self.units:
            units = 'radian'
        else:
            raise ValueError('Could not figure out units {}'.format(
                self.units))
        if self.meshName is None:
            self.meshName = '{}x{}{}'.format(round_res(abs(dLat)),
                                             round_res(abs(dLon)), units)
