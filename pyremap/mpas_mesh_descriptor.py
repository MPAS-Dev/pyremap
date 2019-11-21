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

from pyremap.mesh_descriptor import MeshDescriptor, _create_scrip


class MpasMeshDescriptor(MeshDescriptor):  # {{{
    '''
    A class for describing an MPAS mesh
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, fileName, meshName=None):  # {{{
        '''
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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        ds = xarray.open_dataset(fileName)

        if meshName is None:
            if 'meshName' not in ds.attrs:
                raise ValueError('No meshName provided or found in file.')
            self.meshName = ds.attrs['meshName']
        else:
            self.meshName = meshName

        self.fileName = fileName
        self.regional = True

        # build coords
        self.coords = {'latCell': {'dims': 'nCells',
                                   'data': ds.latCell.values,
                                   'attrs': {'units': 'radians'}},
                       'lonCell': {'dims': 'nCells',
                                   'data': ds.lonCell.values,
                                   'attrs': {'units': 'radians'}}}
        self.dims = ['nCells']
        self.dimSize = [ds.dims[dim] for dim in self.dims]
        ds.close()  # }}}

    def to_scrip(self, scripFileName):  # {{{
        '''
        Given an MPAS mesh file, create a SCRIP file based on the mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.scripFileName = scripFileName

        inFile = netCDF4.Dataset(self.fileName, 'r')
        outFile = netCDF4.Dataset(scripFileName, 'w')

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

        _create_scrip(outFile, grid_size=nCells, grid_corners=maxVertices,
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
        grid_corner_lon = numpy.zeros((nCells, maxVertices))
        grid_corner_lat = numpy.zeros((nCells, maxVertices))
        for iVertex in range(maxVertices):
            cellIndices = numpy.arange(nCells)
            # repeat the last vertex wherever iVertex > nEdgesOnCell
            localVertexIndices = numpy.minimum(nEdgesOnCell - 1, iVertex)
            vertexIndices = verticesOnCell[cellIndices, localVertexIndices] - 1
            grid_corner_lat[cellIndices, iVertex] = latVertex[vertexIndices]
            grid_corner_lon[cellIndices, iVertex] = lonVertex[vertexIndices]

        outFile.variables['grid_corner_lat'][:] = grid_corner_lat[:]
        outFile.variables['grid_corner_lon'][:] = grid_corner_lon[:]

        # Update history attribute of netCDF file
        if hasattr(inFile, 'history'):
            newhist = '\n'.join([getattr(inFile, 'history'),
                                 ' '.join(sys.argv[:])])
        else:
            newhist = sys.argv[:]
        setattr(outFile, 'history', newhist)

        inFile.close()
        outFile.close()  # }}}
    # }}}
