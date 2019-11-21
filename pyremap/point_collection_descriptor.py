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

from pyremap.mesh_descriptor import MeshDescriptor, _create_scrip


class PointCollectionDescriptor(MeshDescriptor):  # {{{
    """
    A class for describing a collection of points

    Author
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/19/2017
    """

    def __init__(self, lats, lons, collectionName,
                 units='degrees', outDimension='nPoints'):  # {{{
        """
        Constructor stores

        Parameters
        ----------

        lats, lons : 1D numpy arrays
            The latitude and longitude of each point

        collectionName : str
            A unique name for the collection of transects, used in the names
            of files containing data mapped to these points.

        units : {'degrees', 'radians'}, optional
            The units of ``lats`` and ``lons``

        outDimension : str, optional
            The name of the dimension corresponding to the points (i.e. the
            "horizontal" dimension of the point collection)

        Author
        ------
        Xylar Asay-Davis

        Last Modified
        -------------
        04/19/2017
        """

        self.meshName = collectionName

        self.regional = True

        self.lat = lats
        self.lon = lons
        self.units = units

        # build coords
        self.coords = {'lat': {'dims': outDimension,
                               'data': self.lat,
                               'attrs': {'units': units}},
                       'lon': {'dims': outDimension,
                               'data': self.lon,
                               'attrs': {'units': units}}}
        self.dims = [outDimension]
        self.dimSize = [len(self.lat)]
        # }}}

    def to_scrip(self, scripFileName):  # {{{
        """
        Given an MPAS mesh file, create a SCRIP file based on the mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        """
        # Authors
        # ------
        # Xylar Asay-Davis

        self.scripFileName = scripFileName

        outFile = netCDF4.Dataset(scripFileName, 'w')

        nPoints = len(self.lat)

        _create_scrip(outFile, grid_size=nPoints,
                      grid_corners=4,
                      grid_rank=1, units=self.units, meshName=self.meshName)

        grid_area = outFile.createVariable('grid_area', 'f8', ('grid_size',))
        grid_area.units = 'radian^2'
        # SCRIP uses square radians
        grid_area[:] = numpy.zeros(nPoints)

        outFile.variables['grid_center_lat'][:] = self.lat
        outFile.variables['grid_center_lon'][:] = self.lon
        outFile.variables['grid_dims'][:] = nPoints
        outFile.variables['grid_imask'][:] = 1

        # grid corners:
        grid_corner_lon = numpy.zeros((nPoints, 4))
        grid_corner_lat = numpy.zeros((nPoints, 4))
        # just repeate the center lat and lon
        for iVertex in range(4):
            grid_corner_lat[:, iVertex] = self.lat
            grid_corner_lon[:, iVertex] = self.lon

        outFile.variables['grid_corner_lat'][:] = grid_corner_lat[:]
        outFile.variables['grid_corner_lon'][:] = grid_corner_lon[:]

        # Update history attribute of netCDF file
        setattr(outFile, 'history', ' '.join(sys.argv[:]))

        outFile.close()  # }}}
# }}}
