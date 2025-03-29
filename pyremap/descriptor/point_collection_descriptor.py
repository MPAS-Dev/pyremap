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
import numpy as np

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import add_history, create_scrip


class PointCollectionDescriptor(MeshDescriptor):
    """
    A class for describing a collection of points

    Attributes
    ----------
    lat : numpy.ndarray
        The latitude of each point

    lon : numpy.ndarray
        The longitude of each point

    units : {'degrees', 'radians'}
        The units of ``lats`` and ``lons``

    history : str
        The history attribute written to SCRIP files
    """

    def __init__(self, lats, lons, collectionName, units='degrees',
                 outDimension='nPoints'):
        """
        Constructor stores

        Parameters
        ----------
        lats : numpy.ndarray
            The latitude of each point

        lons : numpy.ndarray
            The longitude of each point

        collectionName : str
            A unique name for the collection of transects, used in the names
            of files containing data mapped to these points.

        units : {'degrees', 'radians'}, optional
            The units of ``lats`` and ``lons``

        outDimension : str, optional
            The name of the dimension corresponding to the points (i.e. the
            "horizontal" dimension of the point collection)
        """
        super().__init__(meshName=collectionName, regional=True)

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
        self.history = add_history()

    def to_scrip(self, scripFileName, expandDist=None, expandFactor=None):
        """
        Given an MPAS mesh file, create a SCRIP file based on the mesh.

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

        outFile = netCDF4.Dataset(scripFileName, 'w', format=self.format)

        nPoints = len(self.lat)

        create_scrip(outFile, grid_size=nPoints, grid_corners=4,
                     grid_rank=1, units=self.units, meshName=self.meshName)

        grid_area = outFile.createVariable('grid_area', 'f8', ('grid_size',))
        grid_area.units = 'radian^2'
        # SCRIP uses square radians
        grid_area[:] = np.zeros(nPoints)

        outFile.variables['grid_center_lat'][:] = self.lat
        outFile.variables['grid_center_lon'][:] = self.lon
        outFile.variables['grid_dims'][:] = nPoints
        outFile.variables['grid_imask'][:] = 1

        # grid corners:
        grid_corner_lon = np.zeros((nPoints, 4))
        grid_corner_lat = np.zeros((nPoints, 4))
        # just repeat the center lat and lon
        for iVertex in range(4):
            grid_corner_lat[:, iVertex] = self.lat
            grid_corner_lon[:, iVertex] = self.lon

        outFile.variables['grid_corner_lat'][:] = grid_corner_lat[:]
        outFile.variables['grid_corner_lon'][:] = grid_corner_lon[:]

        # Update history attribute of netCDF file
        setattr(outFile, 'history', self.history)

        outFile.close()
