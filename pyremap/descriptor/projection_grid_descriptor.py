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
import pyproj
import xarray as xr

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import (
    add_history,
    expand_scrip,
    interp_extrap_corner,
    unwrap_corners,
    write_netcdf,
)


class ProjectionGridDescriptor(MeshDescriptor):
    """
    A class for describing a general logically rectangular grid that can be
    defined by a ``pyproj`` projection.

    Attributes
    ----------
    projection : pyproj.Proj
        The projection used to map from grid x-y space to latitude and
        longitude

    latLonProjection : pyproj.Proj
        lat-lon projection used to transform from x-y to lat-lon space

    x : numpy.ndarray
        The latitude coordinate at grid-cell centers

    y : numpy.ndarray
        The longitude coordinate at grid-cell centers

    xCorner : numpy.ndarray
        The latitude coordinate at grid-cell corners

    yCorner : numpy.ndarray
        The longitude coordinate at grid-cell corners

    history : str
        The history attribute written to SCRIP files

    xVarName : str
        The name of the x variable

    yVarName : str
        The name of the y variable

    history : str
        The history attribute written to SCRIP files
    """
    def __init__(self, projection, meshName=None):
        """
        Constructor stores the projection

        Parameters
        ----------
        projection : pyproj.Proj
            The projection used to map from grid x-y space to latitude and
            longitude

        meshName : str, optional
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``)
        """
        super().__init__(meshName=meshName, regional=True)

        self.projection = projection
        self.latLonProjection = pyproj.Proj(proj='latlong', datum='WGS84')

        self.x = None
        self.y = None
        self.xCorner = None
        self.yCorner = None
        self.history = None
        self.xVarName = None
        self.yVarName = None
        self.history = None

    @classmethod
    def read(cls, projection, fileName, meshName=None, xVarName='x',
             yVarName='y'):
        """
        Given a grid file with x and y coordinates defining the axes of the
        logically rectangular grid, read in the x and y coordinates and
        interpolate/extrapolate to locate corners.

        Parameters
        ----------
        projection : pyproj.Proj object
            The projection used to map from grid x-y space to latitude and
            longitude

        fileName : str
            The path of the file containing the grid data

        meshName : str, optional
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``).  If not
            provided, the data set in ``fileName`` must have a global
            attribute ``meshName`` that will be used instead.

        xVarName, yVarName : str, optional
            The name of the x and y (in meters) variables in the grid file
        """
        ds = xr.open_dataset(fileName)

        if meshName is None:
            if 'meshName' not in ds.attrs:
                raise ValueError('No meshName provided or found in file.')
            meshName = ds.attrs['meshName']

        descriptor = cls(projection, meshName=meshName)

        # Get info from input file
        descriptor.x = np.array(ds[xVarName].values, float)
        descriptor.y = np.array(ds[yVarName].values, float)

        descriptor._set_coords(xVarName, yVarName, ds[xVarName].dims[0],
                               ds[yVarName].dims[0])

        # interp/extrap corners
        descriptor.xCorner = interp_extrap_corner(descriptor.x)
        descriptor.yCorner = interp_extrap_corner(descriptor.y)

        descriptor.history = add_history(ds=ds)
        return descriptor

    @classmethod
    def create(cls, projection, x, y, meshName):
        """
        Given x and y coordinates defining the axes of the logically
        rectangular grid, save the coordinates interpolate/extrapolate to
        locate corners.

        Parameters
        ----------
        projection : pyproj.Proj object
            The projection used to map from grid x-y space to latitude and
            longitude

        x : numpy.ndarray
            One dimensional array defining the x coordinate of grid cell
            centers.

        y : numpy.ndarray
            One dimensional array defining the y coordinate of grid cell
            centers.

        meshName : str
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``)
        """
        descriptor = cls(projection, meshName=meshName)

        descriptor.x = x
        descriptor.y = y

        descriptor._set_coords('x', 'y', 'x', 'y')

        # interp/extrap corners
        descriptor.xCorner = interp_extrap_corner(descriptor.x)
        descriptor.yCorner = interp_extrap_corner(descriptor.y)
        descriptor.history = add_history()
        return descriptor

    def to_scrip(self, scripFileName, expandDist=None, expandFactor=None):
        """
        Create a SCRIP file based on the grid and projection.

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

        nx = len(self.x)
        ny = len(self.y)

        grid_size = nx * ny

        (center_x, center_y) = np.meshgrid(self.x, self.y)
        (corner_x, corner_y) = np.meshgrid(self.xCorner, self.yCorner)

        (center_lat, center_lon) = self.project_to_lat_lon(center_x, center_y)
        (corner_lat, corner_lon) = self.project_to_lat_lon(corner_x, corner_y)

        ds['grid_center_lat'] = (('grid_size',), center_lat.flat)
        ds['grid_center_lon'] = (('grid_size',), center_lon.flat)
        ds['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'), unwrap_corners(corner_lat)
        )
        ds['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'), unwrap_corners(corner_lon)
        )

        ds['grid_dims'] = xr.DataArray(
            [nx, ny],
            dims=('grid_rank',)
        ).astype('int32')

        ds['grid_imask'] = xr.DataArray(
            np.ones(grid_size, dtype='int32'),
            dims=('grid_size',)
        )

        if expandDist is not None or expandFactor is not None:
            expand_scrip(ds, expandDist, expandFactor)

        ds.grid_center_lat.attrs['units'] = 'degrees'
        ds.grid_center_lon.attrs['units'] = 'degrees'
        ds.grid_corner_lat.attrs['units'] = 'degrees'
        ds.grid_corner_lon.attrs['units'] = 'degrees'
        ds.grid_imask.attrs['units'] = 'unitless'

        ds.attrs['meshName'] = self.meshName
        ds.attrs['history'] = self.history
        write_netcdf(ds, scripFileName, format=self.format)

    def project_to_lat_lon(self, X, Y):
        """
        Given X and Y locations of points in a projection, returns the
        corresponding latitude and longitude of each point.

        Parameters
        ----------
        X : numpy.ndarray
            x array of points in the projection

        Y : numpy.ndarray
            y array of points in the projection

        Returns
        -------
        Lat : numpy.ndarray
            The latitude in degrees with the same size as X and Y

        Lon : numpy.ndarray
            The longitude in degrees with the same size as X and Y
        """
        transformer = pyproj.Transformer.from_proj(self.projection,
                                                   self.latLonProjection)
        Lon, Lat = transformer.transform(X, Y)

        return Lat, Lon

    def _set_coords(self, xVarName, yVarName, xDimName, yDimName):
        """
        Set up a coords dict with x, y, lat and lon
        """
        self.xVarName = xVarName
        self.yVarName = yVarName
        (X, Y) = np.meshgrid(self.x, self.y)
        (Lat, Lon) = self.project_to_lat_lon(X, Y)

        self.coords = {xVarName: {'dims': xDimName,
                                  'data': self.x,
                                  'attrs': {'units': 'meters'}},
                       yVarName: {'dims': yDimName,
                                  'data': self.y,
                                  'attrs': {'units': 'meters'}},
                       'lat': {'dims': (yDimName, xDimName),
                               'data': Lat,
                               'attrs': {'units': 'degrees'}},
                       'lon': {'dims': (yDimName, xDimName),
                               'data': Lon,
                               'attrs': {'units': 'degrees'}}}

        self.dims = [yDimName, xDimName]
        self.dimSize = [len(self.y), len(self.x)]
