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

import numpy
from collections import OrderedDict

from pyremap.mesh_descriptor import MeshDescriptor


class PointCollectionDescriptor(MeshDescriptor):
    """
    A class for describing a collection of points
    """

    def __init__(self, lats, lons, collection_name,
                 units='degrees', out_dimension='nPoints'):
        """
        Constructor stores

        Parameters
        ----------

        lats, lons : numpy array
            The latitude and longitude of each point

        collection_name : str
            A unique name for the collection of transects, used in the names
            of files containing data mapped to these points.

        units : {'degrees', 'radians'}, optional
            The units of ``lats`` and ``lons``

        out_dimension : str, optional
            The name of the dimension corresponding to the points (i.e. the
            "horizontal" dimension of the point collection)
        """

        super().__init__()

        self.meshname = collection_name

        self.regional = True
        # build coords
        self.coords = OrderedDict(
            [('lat', {'dims': out_dimension,
                      'data': lats,
                      'attrs': {'units': units}}),
             ('lon', {'dims': out_dimension,
                      'data': lons,
                      'attrs': {'units': units}})])
        self.lon_lat_coords = ['lat', 'lon']
        self.sizes = OrderedDict([(out_dimension, len(lats))])

        degrees = 'degree' in units
        self.set_lon_lat_cell_centers(lats, lons, degrees=degrees)

        point_count = len(lats)
        lonvert = numpy.zeros((point_count, 4))
        latvert = numpy.zeros((point_count, 4))
        # just repeat the center lat and lon
        for vert_index in range(4):
            lonvert[:, vert_index] = lons
            latvert[:, vert_index] = lats

        vertices_on_cell = numpy.arange(4*point_count).reshape(point_count, 4)

        self.set_lon_lat_vertices(latvert, lonvert,
                                  vertices_on_cell=vertices_on_cell,
                                  degrees=degrees)
