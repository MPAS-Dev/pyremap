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

import netCDF4
import numpy as np
import xarray as xr

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import create_scrip


class MpasEdgeMeshDescriptor(MeshDescriptor):
    """
    A class for describing an MPAS edge mesh

    Attributes
    ----------
    fileName : str
        The path of the file containing the MPAS mesh
    """
    def __init__(self, fileName, meshName=None):
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
        """
        super().__init__()

        with xr.open_dataset(fileName) as ds:

            if meshName is None:
                if 'meshName' not in ds.attrs:
                    raise ValueError('No meshName provided or found in file.')
                self.meshName = ds.attrs['meshName']
            else:
                self.meshName = meshName

            self.fileName = fileName
            self.regional = True

            # build coords
            self.coords = {'latEdge': {'dims': 'nEdges',
                                       'data': ds.latEdge.values,
                                       'attrs': {'units': 'radians'}},
                           'lonEdge': {'dims': 'nEdges',
                                       'data': ds.lonEdge.values,
                                       'attrs': {'units': 'radians'}}}
            self.dims = ['nEdges']
            self.dimSize = [ds.dims[dim] for dim in self.dims]

    def to_scrip(self, scripFileName):
        """
        Given an MPAS mesh file, create a SCRIP file based on the mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        """

        inFile = netCDF4.Dataset(self.fileName, 'r')
        outFile = netCDF4.Dataset(scripFileName, 'w')

        # Get info from input file
        latCell = inFile.variables['latCell'][:]
        lonCell = inFile.variables['lonCell'][:]
        latVertex = inFile.variables['latVertex'][:]
        lonVertex = inFile.variables['lonVertex'][:]
        latEdge = inFile.variables['latEdge'][:]
        lonEdge = inFile.variables['lonEdge'][:]
        verticesOnEdge = inFile.variables['verticesOnEdge'][:]
        cellsOnEdge = inFile.variables['cellsOnEdge'][:]
        nEdges = len(inFile.dimensions['nEdges'])
        dcEdge = inFile.variables['dcEdge'][:]
        dvEdge = inFile.variables['dvEdge'][:]
        sphereRadius = float(inFile.sphere_radius)

        create_scrip(outFile, grid_size=nEdges, grid_corners=4,
                     grid_rank=1, units='radians', meshName=self.meshName)

        validCellsOnEdge = np.zeros(nEdges)
        for iCell in range(2):
            mask = cellsOnEdge[:, iCell] > 0
            validCellsOnEdge[mask] = validCellsOnEdge[mask] + 1

        grid_area = outFile.createVariable('grid_area', 'f8', ('grid_size',))
        grid_area.units = 'radian^2'
        # a triangle if there is 1 valid cell and a diamond if there are 2
        areaEdge = 0.5 * validCellsOnEdge * dcEdge * dvEdge
        # SCRIP uses square radians
        grid_area[:] = areaEdge[:] / (sphereRadius**2)

        outFile.variables['grid_center_lat'][:] = latEdge[:]
        outFile.variables['grid_center_lon'][:] = lonEdge[:]
        outFile.variables['grid_dims'][:] = nEdges
        outFile.variables['grid_imask'][:] = 1

        # grid corners:
        grid_corner_lon = np.zeros((nEdges, 4))
        grid_corner_lat = np.zeros((nEdges, 4))

        # start by repeating vertices, since these always exist
        vertices = verticesOnEdge[:, 0] - 1
        grid_corner_lon[:, 0] = lonVertex[vertices]
        grid_corner_lat[:, 0] = latVertex[vertices]
        grid_corner_lon[:, 1] = lonVertex[vertices]
        grid_corner_lat[:, 1] = latVertex[vertices]
        vertices = verticesOnEdge[:, 1] - 1
        grid_corner_lon[:, 2] = lonVertex[vertices]
        grid_corner_lat[:, 2] = latVertex[vertices]
        grid_corner_lon[:, 3] = lonVertex[vertices]
        grid_corner_lat[:, 3] = latVertex[vertices]

        # replace corners 0 and 2 with cells where they exist
        mask = cellsOnEdge[:, 0] > 0
        cells = cellsOnEdge[mask, 0] - 1
        grid_corner_lon[mask, 0] = lonCell[cells]
        grid_corner_lat[mask, 0] = latCell[cells]
        mask = cellsOnEdge[:, 1] > 0
        cells = cellsOnEdge[mask, 1] - 1
        grid_corner_lon[mask, 2] = lonCell[cells]
        grid_corner_lat[mask, 2] = latCell[cells]

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
        with xr.open_dataset(self.fileName) as ds:

            nCells = ds.sizes['nCells']
            nodeCount = nCells + ds.sizes['nVertices']
            elementCount = ds.sizes['nEdges']
            coordDim = 2

            nodeCoords = np.zeros((nodeCount, coordDim))
            nodeCoords[0:nCells, 0] = ds.lonCell.values
            nodeCoords[0:nCells, 1] = ds.latCell.values
            nodeCoords[nCells:, 0] = ds.lonVertex.values
            nodeCoords[nCells:, 1] = ds.latVertex.values
            centerCoords = np.zeros((elementCount, coordDim))
            centerCoords[:, 0] = ds.lonEdge.values
            centerCoords[:, 1] = ds.latEdge.values

            elementConn = -1 * np.ones((elementCount, 4), dtype=int)

            verticesOnEdge = ds.verticesOnEdge.values
            cellsOnEdge = ds.cellsOnEdge.values

        elementConn[:, 0] = cellsOnEdge[:, 0]
        elementConn[:, 1] = nCells + verticesOnEdge[:, 0]
        elementConn[:, 2] = cellsOnEdge[:, 1]
        elementConn[:, 3] = nCells + verticesOnEdge[:, 1]

        # invalid value is -1, not 0
        elementConn[elementConn == 0] = -1
        for vertex in range(3):
            # swap -1 value until they're at the end
            mask = elementConn[:, vertex] == -1
            elementConn[mask, vertex] = elementConn[mask, vertex + 1]
            elementConn[mask, vertex + 1] = -1

        numElementConn = np.sum((elementConn != -1), axis=1).astype(np.byte)

        ds_out = xr.Dataset()
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
            (('elementCount',), numElementConn)
        ds_out.attrs['gridType'] = 'unstructured mesh'
        ds_out.attrs['version'] = '0.9'

        ds_out.to_netcdf(esmfFileName)
