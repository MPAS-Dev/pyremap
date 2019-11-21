# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

import netCDF4
import numpy
import sys
import pyproj
import xarray

from pyremap.mesh_descriptor import MeshDescriptor, _create_scrip, \
    _unwrap_corners, interp_extrap_corner


class ProjectionGridDescriptor(MeshDescriptor):  # {{{
    '''
    A class for describing a general logically rectangular grid that can be
    defined by a `pyproj` projection.
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, projection):  # {{{
        '''
        Constructor stores the projection

        Parameters
        ----------
        projection : ``pyproj.Proj`` object
            The projection used to map from grid x-y space to latitude and
            longitude
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.projection = projection
        self.latLonProjection = pyproj.Proj(proj='latlong', datum='WGS84')
        self.regional = True

    @classmethod
    def read(cls, projection, fileName, meshName=None, xVarName='x',
             yVarName='y'):  # {{{
        '''
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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        descriptor = cls(projection)

        ds = xarray.open_dataset(fileName)

        if meshName is None:
            if 'meshName' not in ds.attrs:
                raise ValueError('No meshName provided or found in file.')
            descriptor.meshName = ds.attrs['meshName']
        else:
            descriptor.meshName = meshName

        # Get info from input file
        descriptor.x = numpy.array(ds[xVarName].values, float)
        descriptor.y = numpy.array(ds[yVarName].values, float)

        descriptor._set_coords(xVarName, yVarName, ds[xVarName].dims[0],
                               ds[yVarName].dims[0])

        # interp/extrap corners
        descriptor.xCorner = interp_extrap_corner(descriptor.x)
        descriptor.yCorner = interp_extrap_corner(descriptor.y)

        # Update history attribute of netCDF file
        if 'history' in ds.attrs:
            descriptor.history = '\n'.join([ds.attrs['history'],
                                            ' '.join(sys.argv[:])])
        else:
            descriptor.history = sys.argv[:]
        return descriptor  # }}}

    @classmethod
    def create(cls, projection, x, y, meshName):  # {{{
        """
        Given x and y coordinates defining the axes of the logically
        rectangular grid, save the coordinates interpolate/extrapolate to
        locate corners.

        Parameters
        ----------
        projection: pyproj.Proj
            The projection for the grid

        x, y : 1D numpy.arrays
            One dimensional arrays defining the x and y coordinates of grid
            cell centers.

        meshName : str
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        descriptor = cls(projection)

        descriptor.meshName = meshName

        descriptor.x = x
        descriptor.y = y

        descriptor._set_coords('x', 'y', 'x', 'y')

        # interp/extrap corners
        descriptor.xCorner = interp_extrap_corner(descriptor.x)
        descriptor.yCorner = interp_extrap_corner(descriptor.y)
        descriptor.history = sys.argv[:]
        return descriptor  # }}}

    def to_scrip(self, scripFileName):  # {{{
        '''
        Create a SCRIP file based on the grid and projection.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.scripFileName = scripFileName

        outFile = netCDF4.Dataset(scripFileName, 'w')

        nx = len(self.x)
        ny = len(self.y)

        grid_size = nx * ny

        _create_scrip(outFile, grid_size=grid_size, grid_corners=4,
                      grid_rank=2, units='degrees', meshName=self.meshName)

        (X, Y) = numpy.meshgrid(self.x, self.y)
        (XCorner, YCorner) = numpy.meshgrid(self.xCorner, self.yCorner)

        (Lat, Lon) = self.project_to_lat_lon(X, Y)
        (LatCorner, LonCorner) = self.project_to_lat_lon(XCorner, YCorner)

        outFile.variables['grid_center_lat'][:] = Lat.flat
        outFile.variables['grid_center_lon'][:] = Lon.flat
        outFile.variables['grid_dims'][:] = [nx, ny]
        outFile.variables['grid_imask'][:] = 1

        outFile.variables['grid_corner_lat'][:] = _unwrap_corners(LatCorner)
        outFile.variables['grid_corner_lon'][:] = _unwrap_corners(LonCorner)

        setattr(outFile, 'history', self.history)

        outFile.close()  # }}}

    def project_to_lat_lon(self, X, Y):  # {{{
        '''
        Given X and Y locations of points in a projection, returns the
        corresponding latitude and longitude of each point.

        Parameters
        ----------
        outFile : file pointer
            A SCRIP file opened in write mode

        X, Y : 1D or 2D numpy.array
            X and y arrays of points in in the projection

        Returns
        -------
        Lat, Lon : numpy.array with same shape as X and Y
            the latitude and longitude in degrees of the points
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        Lon, Lat = pyproj.transform(self.projection, self.latLonProjection,
                                    X, Y)

        return (Lat, Lon)  # }}}

    def _set_coords(self, xVarName, yVarName, xDimName, yDimName):  # {{{
        '''
        Set up a coords dict with x, y, lat and lon
        '''
        self.xVarName = xVarName
        self.yVarName = yVarName
        (X, Y) = numpy.meshgrid(self.x, self.y)
        (Lat, Lon) = self.project_to_lat_lon(X, Y)

        self.coords = {xVarName: {'dims': xDimName,
                                  'data': self.x,
                                  'attrs': {'units': 'meters'}},
                       yVarName: {'dims': yDimName,
                                  'data': self.y,
                                  'attrs': {'units': 'meters'}},
                       'lat': {'dims': (xDimName, yDimName),
                               'data': Lat,
                               'attrs': {'units': 'degrees'}},
                       'lon': {'dims': (xDimName, yDimName),
                               'data': Lon,
                               'attrs': {'units': 'degrees'}}}

        self.dims = [xDimName, yDimName]
        self.dimSize = [len(self.x), len(self.y)]
        # }}}

    # }}}
