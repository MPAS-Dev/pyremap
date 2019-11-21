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
import xarray

from pyremap.mesh_descriptor import MeshDescriptor, _create_scrip, \
    _unwrap_corners, interp_extrap_corner, _round_res


def get_lat_lon_descriptor(dLon, dLat, lonMin=-180., lonMax=180., latMin=-90.,
                           latMax=90.):
    """
    Get a descriptor of a lat/lon grid, used for remapping

    Parameters
    ----------
    dLon, dLat :  float`
        Resolution of the lon-lat grid in degrees

    Returns
    -------
    descriptor : ``LatLonGridDescriptor`` object
        A descriptor of the lat/lon grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    nLat = int((latMax - latMin) / dLat) + 1
    nLon = int((lonMax - lonMin) / dLon) + 1
    lat = numpy.linspace(latMin, latMax, nLat)
    lon = numpy.linspace(lonMin, lonMax, nLon)

    descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

    return descriptor


class LatLonGridDescriptor(MeshDescriptor):  # {{{
    '''
    A class for describing a lat-lon grid
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self):  # {{{
        '''
        Constructor stores the file name

        Parameters
        ----------
        fileName : str
            The path of the file containing the MPAS mesh
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.regional = False
        self.meshName = None  # }}}

    @classmethod
    def read(cls, fileName=None, ds=None, latVarName='lat',
             lonVarName='lon'):  # {{{
        '''
        Read the lat-lon grid from a file with the given lat/lon var names.

        Parameters
        ----------
        fileName : str, optional
            The path of the file containing the lat-lon grid (if ``ds`` is not
            supplied directly)

        ds : ``xarray.Dataset`` object, optional
            The path of the file containing the lat-lon grid (if supplied,
            ``fileName`` will be ignored)

        latVarName, lonVarName : str, optional
            The name of the latitude and longitude variables in the grid file
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        if ds is None:
            ds = xarray.open_dataset(fileName)

        descriptor = cls()

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
        return descriptor  # }}}

    @classmethod
    def create(cls, latCorner, lonCorner, units='degrees'):  # {{{
        '''
        Create the lat-lon grid with the given arrays and units.

        Parameters
        ----------
        latCorner, lonCorner : 1D numpy.arrays
            One dimensional arrays defining the latitude and longitude
            coordinates of grid corners.

        units : {'degrees', 'radians'}, optional
            The units of `latCorner` and `lonCorner`
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        descriptor = cls()

        descriptor.latCorner = latCorner
        descriptor.lonCorner = lonCorner
        descriptor.lon = 0.5 * (lonCorner[0:-1] + lonCorner[1:])
        descriptor.lat = 0.5 * (latCorner[0:-1] + latCorner[1:])
        descriptor.units = units
        descriptor.history = sys.argv[:]
        descriptor._set_coords('lat', 'lon', 'lat', 'lon')
        return descriptor  # }}}

    def to_scrip(self, scripFileName):  # {{{
        '''
        Given a lat-lon grid file, create a SCRIP file based on the grid.

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

        nLat = len(self.lat)
        nLon = len(self.lon)

        grid_size = nLat * nLon

        _create_scrip(outFile, grid_size=grid_size, grid_corners=4,
                      grid_rank=2, units=self.units, meshName=self.meshName)

        (Lon, Lat) = numpy.meshgrid(self.lon, self.lat)
        (LonCorner, LatCorner) = numpy.meshgrid(self.lonCorner, self.latCorner)

        outFile.variables['grid_center_lat'][:] = Lat.flat
        outFile.variables['grid_center_lon'][:] = Lon.flat
        outFile.variables['grid_dims'][:] = [nLon, nLat]
        outFile.variables['grid_imask'][:] = 1

        outFile.variables['grid_corner_lat'][:] = _unwrap_corners(LatCorner)
        outFile.variables['grid_corner_lon'][:] = _unwrap_corners(LonCorner)

        setattr(outFile, 'history', self.history)

        outFile.close()  # }}}

    def _set_coords(self, latVarName, lonVarName, latDimName,
                    lonDimName):  # {{{
        '''
        Set up a coords dict with lat and lon
        '''
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
            self.meshName = '{}x{}{}'.format(_round_res(abs(dLat)),
                                             _round_res(abs(dLon)), units)

        # }}}
    # }}}
