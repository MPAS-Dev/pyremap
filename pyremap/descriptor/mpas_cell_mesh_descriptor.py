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

import warnings

import netCDF4
import numpy as np
import xarray as xr

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import add_history, create_scrip, expand_scrip


class MpasCellMeshDescriptor(MeshDescriptor):
    """
    A class for describing an MPAS cell mesh

    Attributes
    ----------
    vertices : bool
        Whether the mapping is to or from vertices instead of corners
        (for non-conservative remapping)

    fileName : str
        The path of the file containing the MPAS mesh

    history : str
        The history attribute written to SCRIP files
    """
    def __init__(self, fileName, meshName=None, vertices=False):
        """
        Constructor stores the file name

        Parameters
        ----------
        fileName : str
            The path of the file containing the MPAS mesh

        meshName : str, optional
            The name of the MPAS mesh (e.g. ``'oEC60to30'`` or
            ``'oRRS18to6'``).  If not provided, the data set in ``fileName``
            must have a global attribute ``meshName`` that will be used
            instead.

        vertices : bool, optional
            Whether the mapping is to or from vertices instead of corners
            (for non-conservative remapping)
        """
        super().__init__()

        if vertices:
            warnings.warn('Creating a MpasCellMeshDescriptor with '
                          'vertices=True is deprecated and will be removed in '
                          'the next release. Use MpasVertexMeshDescriptor '
                          'instead.', DeprecationWarning)

        self.vertices = vertices

        with xr.open_dataset(fileName) as ds:

            if meshName is None:
                if 'meshName' not in ds.attrs:
                    raise ValueError('No meshName provided or found in file.')
                self.meshName = ds.attrs['meshName']
            else:
                self.meshName = meshName

            self.fileName = fileName
            self.regional = True

            if vertices:
                # build coords
                self.coords = {'latVertex': {'dims': 'nVertices',
                                             'data': ds.latVertex.values,
                                             'attrs': {'units': 'radians'}},
                               'lonVertex': {'dims': 'nVertices',
                                             'data': ds.lonVertex.values,
                                             'attrs': {'units': 'radians'}}}
                self.dims = ['nVertices']
            else:
                # build coords
                self.coords = {'latCell': {'dims': 'nCells',
                                           'data': ds.latCell.values,
                                           'attrs': {'units': 'radians'}},
                               'lonCell': {'dims': 'nCells',
                                           'data': ds.lonCell.values,
                                           'attrs': {'units': 'radians'}}}
                self.dims = ['nCells']
            self.dimSize = [ds.sizes[dim] for dim in self.dims]

            self.history = add_history(ds=ds)

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
        if self.vertices:
            raise ValueError('A SCRIP file won\'t work for remapping vertices')

        inFile = netCDF4.Dataset(self.fileName, 'r')
        outFile = netCDF4.Dataset(scripFileName, 'w', format=self.format)

        # Get info from input file
        latCell = inFile.variables['latCell'][:]
        lonCell = inFile.variables['lonCell'][:]
        latVertex = inFile.variables['latVertex'][:]
        lonVertex = inFile.variables['lonVertex'][:]
        verticesOnCell = inFile.variables['verticesOnCell'][:]
        nEdgesOnCell = inFile.variables['nEdgesOnCell'][:]
        nCells = len(inFile.dimensions['nCells'])
        maxVertices = len(inFile.dimensions['maxEdges'])
        areaCell = inFile.variables['areaCell'][:]
        sphereRadius = float(inFile.sphere_radius)

        create_scrip(outFile, grid_size=nCells, grid_corners=maxVertices,
                     grid_rank=1, units='radians', meshName=self.meshName)

        grid_area = outFile.createVariable('grid_area', 'f8', ('grid_size',))
        grid_area.units = 'radian^2'
        # SCRIP uses square radians
        grid_area[:] = areaCell[:] / (sphereRadius**2)

        outFile.variables['grid_center_lat'][:] = latCell[:]
        outFile.variables['grid_center_lon'][:] = lonCell[:]
        outFile.variables['grid_dims'][:] = nCells
        outFile.variables['grid_imask'][:] = 1

        # grid corners:
        grid_corner_lon = np.zeros((nCells, maxVertices))
        grid_corner_lat = np.zeros((nCells, maxVertices))
        for iVertex in range(maxVertices):
            cellIndices = np.arange(nCells)
            # repeat the last vertex wherever iVertex > nEdgesOnCell
            localVertexIndices = np.minimum(nEdgesOnCell - 1, iVertex)
            vertexIndices = verticesOnCell[cellIndices, localVertexIndices] - 1
            grid_corner_lat[cellIndices, iVertex] = latVertex[vertexIndices]
            grid_corner_lon[cellIndices, iVertex] = lonVertex[vertexIndices]

        outFile.variables['grid_corner_lat'][:] = grid_corner_lat[:]
        outFile.variables['grid_corner_lon'][:] = grid_corner_lon[:]

        if expandDist is not None or expandFactor is not None:
            expand_scrip(outFile, expandDist, expandFactor)

        setattr(outFile, 'history', self.history)

        inFile.close()
        outFile.close()
