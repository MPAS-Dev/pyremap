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

import xarray
from collections import OrderedDict

from pyremap.mesh_descriptor import MeshDescriptor


class MpasMeshDescriptor(MeshDescriptor):
    """
    A class for describing an MPAS mesh
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, filename, meshname=None):
        """
        Constructor stores the file name

        Parameters
        ----------
        filename : str
            The path of the file containing the MPAS mesh

        meshname : str, optional
            The name of the MPAS mesh (e.g. ``'oEC60to30'`` or
            ``'oRRS18to6'``).  If not provided, the data set in ``filename``
            must have a global attribute ``meshname`` that will be used
            instead.
        """

        super().__init__()

        with xarray.open_dataset(filename) as ds:

            if meshname is None:
                if 'meshname' not in ds.attrs:
                    raise ValueError('No meshname provided or found in file.')
                self.meshname = ds.attrs['meshname']
            else:
                self.meshname = meshname

            self.regional = True

            self.coords = OrderedDict()
            self.coords['latCell'] = {'dims': ('nCells',),
                                      'data': ds.latCell.values,
                                      'attrs': {'units': 'radians'}}
            self.coords['lonCell'] = {'dims': ('nCells',),
                                      'data': ds.lonCell.values,
                                      'attrs': {'units': 'radians'}}
            self.coords['xCell'] = {'dims': ('nCells',),
                                    'data': ds.xCell.values,
                                    'attrs': {'units': 'meters'}}
            self.coords['yCell'] = {'dims': ('nCells',),
                                    'data': ds.yCell.values,
                                    'attrs': {'units': 'meters'}}
            self.coords['zCell'] = {'dims': ('nCells',),
                                    'data': ds.zCell.values,
                                    'attrs': {'units': 'meters'}}

            self.sizes = OrderedDict([('nCells', ds.sizes['nCells'])])

        vertices_on_cell = ds.verticesOnCell.values - 1
        vertex_count_on_cells = ds.nEdgesOnCell.values

        self.set_lon_lat_vertices(
            ds.lonVertex.values, ds.latVertex.values,
            vertices_on_cell=vertices_on_cell,
            vertex_count_on_cells=vertex_count_on_cells,
            degrees=False)

        self.set_cartesian_vertices(
            ds.xVertex.values, ds.yVertex.values, ds.zVertex.values,
            vertices_on_cell=vertices_on_cell,
            vertex_count_on_cells=vertex_count_on_cells)

        self.set_lon_lat_cell_centers(ds.lonCell.values, ds.latCell.values,
                                      degrees=False)

        self.set_cartesian_cell_centers(ds.xCell.values, ds.yCell.values,
                                        ds.zCell.values)
