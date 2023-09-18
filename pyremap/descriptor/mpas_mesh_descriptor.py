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

from pyremap.descriptor.mpas_cell_mesh_descriptor import MpasCellMeshDescriptor


class MpasMeshDescriptor(MpasCellMeshDescriptor):
    """
    A class for describing an MPAS cell mesh (which also supports
    non-conservative remapping from vertices)
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
        warnings.warn('MpasMeshDescriptor is deprecated and will be removed '
                      'in the next release. Use MpasCellMeshDescriptor '
                      'instead.', DeprecationWarning)

        super().__init__(fileName, meshName=meshName, vertices=vertices)
