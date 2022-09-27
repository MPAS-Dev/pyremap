# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

import netCDF4
import numpy
import sys
import xarray

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import interp_extrap_corner, \
    create_scrip, unwrap_corners, round_res


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
    descriptor : xarray.LatLonGridDescriptor object
        A descriptor of the lat/lon grid
    """
    nLat = int((latMax - latMin) / dLat) + 1
    nLon = int((lonMax - lonMin) / dLon) + 1
    lat = numpy.linspace(latMin, latMax, nLat)
    lon = numpy.linspace(lonMin, lonMax, nLon)

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
    def __init__(self, meshName=None, regional=False):
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
             lonVarName='lon', meshName=None, regional=False):
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
            Whether this is a regional or global grid
        """
        if ds is None:
            ds = xarray.open_dataset(fileName)

        descriptor = cls(meshName=meshName, regional=regional)

        if descriptor.meshName is None and 'meshName' in ds.attrs:
            descriptor.meshName = ds.attrs['meshName']

        # Get info from input file
        descriptor.lat = numpy.array(ds[latVarName].values, float)
        descriptor.lon = numpy.array(ds[lonVarName].values, float)
        if 'degree' in ds[latVarName].units:
            descriptor.units = 'degrees'
        else:
            descriptor.units = 'radians'

        # interp/extrap corners
        descriptor.lonCorner = interp_extrap_corner(descriptor.lon)
        descriptor.latCorner = interp_extrap_corner(descriptor.lat)

        descriptor._set_coords(latVarName, lonVarName, ds[latVarName].dims[0],
                               ds[lonVarName].dims[0])

        if 'history' in ds.attrs:
            descriptor.history = '\n'.join([ds.attrs['history'],
                                            ' '.join(sys.argv[:])])
        else:
            descriptor.history = sys.argv[:]
        return descriptor

    @classmethod
    def create(cls, latCorner, lonCorner, units='degrees'):
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
        """
        descriptor = cls()

        descriptor.latCorner = latCorner
        descriptor.lonCorner = lonCorner
        descriptor.lon = 0.5 * (lonCorner[0:-1] + lonCorner[1:])
        descriptor.lat = 0.5 * (latCorner[0:-1] + latCorner[1:])
        descriptor.units = units
        descriptor.history = sys.argv[:]
        descriptor._set_coords('lat', 'lon', 'lat', 'lon')
        return descriptor

    def to_scrip(self, scripFileName):
        """
        Given a lat-lon grid file, create a SCRIP file based on the grid.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        """
        outFile = netCDF4.Dataset(scripFileName, 'w')

        nLat = len(self.lat)
        nLon = len(self.lon)

        grid_size = nLat * nLon

        create_scrip(outFile, grid_size=grid_size, grid_corners=4,
                     grid_rank=2, units=self.units, meshName=self.meshName)

        (Lon, Lat) = numpy.meshgrid(self.lon, self.lat)
        (LonCorner, LatCorner) = numpy.meshgrid(self.lonCorner, self.latCorner)

        outFile.variables['grid_center_lat'][:] = Lat.flat
        outFile.variables['grid_center_lon'][:] = Lon.flat
        outFile.variables['grid_dims'][:] = [nLon, nLat]
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
        nLat = len(self.lat)
        nLon = len(self.lon)

        (Lon, Lat) = numpy.meshgrid(self.lon, self.lat)
        (LonCorner, LatCorner) = numpy.meshgrid(self.lonCorner, self.latCorner)
        (XIndices, YIndices) = numpy.meshgrid(numpy.arange(nLon + 1),
                                              numpy.arange(nLat + 1))

        elementCount = nLon * nLat
        nodeCount = (nLon + 1) * (nLat + 1)
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
                             (nLon + 1) * YIndices[0:-1, 0:-1].ravel())
        elementConn[:, 1] = (XIndices[0:-1, 1:].ravel() +
                             (nLon + 1) * YIndices[0:-1, 1:].ravel())
        elementConn[:, 2] = (XIndices[1:, 1:].ravel() +
                             (nLon + 1) * YIndices[1:, 1:].ravel())
        elementConn[:, 3] = (XIndices[1:, 0:-1].ravel() +
                             (nLon + 1) * YIndices[1:, 0:-1].ravel())

        ds_out = xarray.Dataset()
        ds_out['nodeCoords'] = (('nodeCount', 'coordDim'), nodeCoords)
        ds_out.nodeCoords.attrs['units'] = self.units
        ds_out['centerCoords'] = \
            (('elementCount', 'coordDim'), centerCoords)
        ds_out.centerCoords.attrs['units'] = self.units
        ds_out['elementConn'] = \
            (('elementCount', 'maxNodePElement'), elementConn + 1)
        ds_out.elementConn.attrs['start_index'] = 1
        ds_out.elementConn.attrs['_FillValue'] = -1
        ds_out['numElementConn'] = \
            (('elementCount',), maxNodePElement * numpy.ones(elementCount,
                                                             dtype=int))

        ds_out['origGridDims'] = (('origGridRank',), [nLon, nLat])

        ds_out.attrs['gridType'] = 'unstructured mesh'
        ds_out.attrs['version'] = '0.9'

        ds_out.to_netcdf(esmfFileName)

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
            if numpy.abs(lonRange - 360.) > 1e-10:
                self.regional = True
            if numpy.abs(latRange - 180.) > 1e-10:
                self.regional = True
        elif 'rad' in self.units:
            if numpy.abs(lonRange - 2. * numpy.pi) > 1e-10:
                self.regional = True
            if numpy.abs(latRange - numpy.pi) > 1e-10:
                self.regional = True
            units = 'radian'
        else:
            raise ValueError('Could not figure out units {}'.format(
                self.units))
        if self.meshName is None:
            self.meshName = '{}x{}{}'.format(round_res(abs(dLat)),
                                             round_res(abs(dLon)), units)
