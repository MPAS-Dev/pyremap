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

import netCDF4
import numpy
import pyproj
import xarray

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import (
    add_history,
    create_scrip,
    interp_extrap_corner,
    unwrap_corners,
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
        ds = xarray.open_dataset(fileName)

        if meshName is None:
            if 'meshName' not in ds.attrs:
                raise ValueError('No meshName provided or found in file.')
            meshName = ds.attrs['meshName']

        descriptor = cls(projection, meshName=meshName)

        # Get info from input file
        descriptor.x = numpy.array(ds[xVarName].values, float)
        descriptor.y = numpy.array(ds[yVarName].values, float)

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

    def to_scrip(self, scripFileName):
        """
        Create a SCRIP file based on the grid and projection.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        """
        outFile = netCDF4.Dataset(scripFileName, 'w', format=self.format)

        nx = len(self.x)
        ny = len(self.y)

        grid_size = nx * ny

        create_scrip(outFile, grid_size=grid_size, grid_corners=4,
                     grid_rank=2, units='degrees', meshName=self.meshName)

        (X, Y) = numpy.meshgrid(self.x, self.y)
        (XCorner, YCorner) = numpy.meshgrid(self.xCorner, self.yCorner)

        (Lat, Lon) = self.project_to_lat_lon(X, Y)
        (LatCorner, LonCorner) = self.project_to_lat_lon(XCorner, YCorner)

        outFile.variables['grid_center_lat'][:] = Lat.flat
        outFile.variables['grid_center_lon'][:] = Lon.flat
        outFile.variables['grid_dims'][:] = [nx, ny]
        outFile.variables['grid_imask'][:] = 1

        outFile.variables['grid_corner_lat'][:] = unwrap_corners(LatCorner)
        outFile.variables['grid_corner_lon'][:] = unwrap_corners(LonCorner)

        setattr(outFile, 'history', self.history)

        outFile.close()

    def to_esmf(self, esmfFileName):
        """
        Create an ESMF mesh file for the mesh

        Parameters
        ----------
        esmfFileName : str
            The path to which the ESMF mesh file should be written
        """
        nx = len(self.x)
        ny = len(self.y)

        (X, Y) = numpy.meshgrid(self.x, self.y)
        (XCorner, YCorner) = numpy.meshgrid(self.xCorner, self.yCorner)
        (XIndices, YIndices) = numpy.meshgrid(numpy.arange(nx + 1),
                                              numpy.arange(ny + 1))

        (Lat, Lon) = self.project_to_lat_lon(X, Y)
        (LatCorner, LonCorner) = self.project_to_lat_lon(XCorner, YCorner)

        elementCount = nx * ny
        nodeCount = (nx + 1) * (ny + 1)
        coordDim = 2
        maxNodePElement = 4

        nodeCoords = numpy.zeros((nodeCount, coordDim))
        nodeCoords[:, 0] = LonCorner.flat
        nodeCoords[:, 1] = LatCorner.flat

        centerCoords = numpy.zeros((elementCount, coordDim))
        centerCoords[:, 0] = Lon.flat
        centerCoords[:, 1] = Lat.flat

        elementConn = numpy.zeros((elementCount, maxNodePElement), dtype=int)
        elementConn[:, 0] = (XIndices[0:-1, 0:-1].ravel() +
                             (nx + 1) * YIndices[0:-1, 0:-1].ravel())
        elementConn[:, 1] = (XIndices[0:-1, 1:].ravel() +
                             (nx + 1) * YIndices[0:-1, 1:].ravel())
        elementConn[:, 2] = (XIndices[1:, 1:].ravel() +
                             (nx + 1) * YIndices[1:, 1:].ravel())
        elementConn[:, 3] = (XIndices[1:, 0:-1].ravel() +
                             (nx + 1) * YIndices[1:, 0:-1].ravel())

        ds_out = xarray.Dataset()
        ds_out['nodeCoords'] = (('nodeCount', 'coordDim'), nodeCoords)
        ds_out.nodeCoords.attrs['units'] = 'degrees'
        ds_out['centerCoords'] = \
            (('elementCount', 'coordDim'), centerCoords)
        ds_out.centerCoords.attrs['units'] = 'degrees'
        ds_out['elementConn'] = \
            (('elementCount', 'maxNodePElement'), elementConn + 1)
        ds_out.elementConn.attrs['start_index'] = 1
        ds_out.elementConn.attrs['_FillValue'] = -1
        ds_out['numElementConn'] = \
            (('elementCount',), maxNodePElement * numpy.ones(elementCount,
                                                             dtype=int))

        ds_out['origGridDims'] = (('origGridRank',), [ny, nx])

        ds_out.attrs['gridType'] = 'unstructured mesh'
        ds_out.attrs['version'] = '0.9'

        ds_out.to_netcdf(esmfFileName, format=self.format,
                         engine=self.engine)

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
        (X, Y) = numpy.meshgrid(self.x, self.y)
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
