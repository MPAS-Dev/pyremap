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

from pyremap.mesh_descriptor import MeshDescriptor, interp_extrap_corners_2d, \
    _create_scrip, _unwrap_corners, _round_res


class LatLon2DGridDescriptor(MeshDescriptor):  # {{{
    """
    A class for describing a lat-lon grid that may not be a tensor grid
    (lat/lon are 2D arrays).  The grid is assumed to be regional, since this
    is difficult to determine just from the lat/lon values.  The calling code
    should set ``regional = False`` for global grids with 2D lat/lon
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self):  # {{{
        """
        Constructor stores the file name

        Parameters
        ----------
        fileName : str
            The path of the file containing the MPAS mesh
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.regional = True
        self.meshName = None  # }}}

    @classmethod
    def read(cls, fileName=None, ds=None, latVarName='lat',
             lonVarName='lon'):  # {{{
        """
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
        """
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
        descriptor.lonCorner = interp_extrap_corners_2d(descriptor.lon)
        descriptor.latCorner = interp_extrap_corners_2d(descriptor.lat)

        descriptor._set_coords(latVarName, lonVarName, ds[latVarName].dims[0],
                               ds[latVarName].dims[1])

        if 'history' in ds.attrs:
            descriptor.history = '\n'.join([ds.attrs['history'],
                                            ' '.join(sys.argv[:])])
        else:
            descriptor.history = sys.argv[:]
        return descriptor  # }}}

    def to_scrip(self, scripFileName):  # {{{
        """
        Given a lat-lon grid file, create a SCRIP file based on the grid.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.scripFileName = scripFileName

        outFile = netCDF4.Dataset(scripFileName, 'w')

        nLat, nLon = self.lat.shape

        grid_size = nLat * nLon

        _create_scrip(outFile, grid_size=grid_size, grid_corners=4,
                      grid_rank=2, units=self.units, meshName=self.meshName)

        outFile.variables['grid_center_lat'][:] = self.lat.flat
        outFile.variables['grid_center_lon'][:] = self.lon.flat
        outFile.variables['grid_dims'][:] = [nLon, nLat]
        outFile.variables['grid_imask'][:] = 1

        outFile.variables['grid_corner_lat'][:] = \
            _unwrap_corners(self.lonCorner)
        outFile.variables['grid_corner_lon'][:] = \
            _unwrap_corners(self.latCorner)

        setattr(outFile, 'history', self.history)

        outFile.close()  # }}}

    def _set_coords(self, latVarName, lonVarName, latDimName,
                    lonDimName):  # {{{
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
            self.meshName = '{}x{}{}'.format(_round_res(abs(dLat)),
                                             _round_res(abs(dLon)), units)

        # }}}
    # }}}


# vim: ai ts=4 sts=4 et sw=4 ft=python
