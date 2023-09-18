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
import xarray as xr

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import add_history, create_scrip


class MpasVertexMeshDescriptor(MeshDescriptor):
    """
    A class for describing an MPAS vertex mesh

    Attributes
    ----------
    fileName : str
        The path of the file containing the MPAS mesh

    history : str
        The history attribute written to SCRIP files
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
            self.coords = {'latVertex': {'dims': 'nVertices',
                                         'data': ds.latVertex.values,
                                         'attrs': {'units': 'radians'}},
                           'lonVertex': {'dims': 'nVertices',
                                         'data': ds.lonVertex.values,
                                         'attrs': {'units': 'radians'}}}
            self.dims = ['nVertices']
            self.dimSize = [ds.dims[dim] for dim in self.dims]

            self.history = add_history(ds=ds)

    def to_scrip(self, scripFileName):
        """
        Given an MPAS mesh file, create a SCRIP file based on the mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        """

        inFile = netCDF4.Dataset(self.fileName, 'r')
        outFile = netCDF4.Dataset(scripFileName, 'w', format=self.format)

        # Get info from input file
        latCell = inFile.variables['latCell'][:]
        lonCell = inFile.variables['lonCell'][:]
        latVertex = inFile.variables['latVertex'][:]
        lonVertex = inFile.variables['lonVertex'][:]
        latEdge = inFile.variables['latEdge'][:]
        lonEdge = inFile.variables['lonEdge'][:]
        cellsOnVertex = inFile.variables['cellsOnVertex'][:]
        edgesOnVertex = inFile.variables['edgesOnVertex'][:]
        nVertices = len(inFile.dimensions['nVertices'])
        vertexDegree = len(inFile.dimensions['vertexDegree'])
        kiteAreasOnVertex = inFile.variables['kiteAreasOnVertex'][:]
        sphereRadius = float(inFile.sphere_radius)

        if vertexDegree != 3:
            raise ValueError(f'MpasVertexMeshDescriptor does not support '
                             f'vertexDegree {vertexDegree}')

        create_scrip(outFile, grid_size=nVertices, grid_corners=6,
                     grid_rank=1, units='radians', meshName=self.meshName)

        validCellsOnVertex = np.zeros(nVertices, dtype=int)
        vertexArea = np.zeros(nVertices)
        for iCell in range(vertexDegree):
            mask = cellsOnVertex[:, iCell] > 0
            validCellsOnVertex[mask] = validCellsOnVertex[mask] + 1
            vertexArea[mask] = (vertexArea[mask] +
                                kiteAreasOnVertex[mask, iCell])

        grid_area = outFile.createVariable('grid_area', 'f8', ('grid_size',))
        grid_area.units = 'radian^2'
        # SCRIP uses square radians
        grid_area[:] = vertexArea[:] / (sphereRadius**2)

        outFile.variables['grid_center_lat'][:] = latVertex[:]
        outFile.variables['grid_center_lon'][:] = lonVertex[:]
        outFile.variables['grid_dims'][:] = nVertices
        outFile.variables['grid_imask'][:] = 1

        # grid corners:
        grid_corner_lon = np.zeros((nVertices, 6))
        grid_corner_lat = np.zeros((nVertices, 6))

        # start by filling the corners with the vertex lat and lon as the
        # fallback
        for iCorner in range(6):
            grid_corner_lon[:, iCorner] = lonVertex
            grid_corner_lat[:, iCorner] = latVertex

        # if edges on vertex exist, replace even indices with their locations
        for iEdge in range(3):
            mask = edgesOnVertex[:, iEdge] > 0
            edges = edgesOnVertex[mask, iEdge] - 1
            grid_corner_lon[mask, 2 * iEdge] = lonEdge[edges]
            grid_corner_lat[mask, 2 * iEdge] = latEdge[edges]

        # similarly, if cells on vertex exist, replace odd indices with their
        # locations
        for iCell in range(3):
            mask = cellsOnVertex[:, iCell] > 0
            cells = cellsOnVertex[mask, iCell] - 1
            grid_corner_lon[mask, 2 * iCell + 1] = lonCell[cells]
            grid_corner_lat[mask, 2 * iCell + 1] = latCell[cells]

        outFile.variables['grid_corner_lat'][:] = grid_corner_lat[:]
        outFile.variables['grid_corner_lon'][:] = grid_corner_lon[:]

        setattr(outFile, 'history', self.history)

        inFile.close()
        outFile.close()
