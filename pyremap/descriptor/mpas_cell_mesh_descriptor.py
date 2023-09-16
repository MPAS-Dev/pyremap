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

import sys
import warnings

import netCDF4
import numpy
import xarray

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import create_scrip


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

        with xarray.open_dataset(fileName) as ds:

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
            self.dimSize = [ds.dims[dim] for dim in self.dims]

    def to_scrip(self, scripFileName):
        """
        Given an MPAS mesh file, create a SCRIP file based on the mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
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
        outFile.close()

    def to_esmf(self, esmfFileName):
        """
        Create an ESMF mesh file for the mesh

        Parameters
        ----------
        esmfFileName : str
            The path to which the ESMF mesh file should be written
        """
        with xarray.open_dataset(self.fileName) as ds:

            nodeCount = ds.sizes['nVertices']
            elementCount = ds.sizes['nCells']
            coordDim = 2

            nodeCoords = numpy.zeros((nodeCount, coordDim))
            nodeCoords[:, 0] = ds.lonVertex.values
            nodeCoords[:, 1] = ds.latVertex.values
            centerCoords = numpy.zeros((elementCount, coordDim))
            centerCoords[:, 0] = ds.lonCell.values
            centerCoords[:, 1] = ds.latCell.values

            elementConn = ds.verticesOnCell.values
            elementConn[elementConn == nodeCount + 1] = -1

            ds_out = xarray.Dataset()
            ds_out['nodeCoords'] = (('nodeCount', 'coordDim'), nodeCoords)
            ds_out.nodeCoords.attrs['units'] = 'radians'
            ds_out['centerCoords'] = \
                (('elementCount', 'coordDim'), centerCoords)
            ds_out.centerCoords.attrs['units'] = 'radians'
            ds_out['elementConn'] = \
                (('elementCount', 'maxNodePElement'), elementConn)
            ds_out.elementConn.attrs['start_index'] = 1
            ds_out.elementConn.attrs['_FillValue'] = -1
            ds_out['numElementConn'] = \
                (('elementCount',), ds.nEdgesOnCell.values)
            ds_out.attrs['gridType'] = 'unstructured mesh'
            ds_out.attrs['version'] = '0.9'

            ds_out.to_netcdf(esmfFileName, format=self.format,
                             engine=self.engine)
