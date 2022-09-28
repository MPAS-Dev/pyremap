# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE


class MeshDescriptor:
    """
    A class for describing a mesh

    Attributes
    ----------
    meshName : str
        The name of the mesh or grid, used to give mapping files unique names.

    regional : bool
        Whether this is a regional or global grid

    dims : list of str
        Dimension names in the mesh or grid

    dimSizes : list of int
        Size of each dimension in ``dims``

    coords : dict
        A dictionary that can be used to construct an xarray.DataArray for each
        coordinate in the mesh or grid
    """
    def __init__(self, meshName=None, regional=None):
        """
        Construct a mesh descriptor

        meshName : str or None, optional
            The name of the mesh or grid, used to give mapping files unique
            names.  If not provided here, it should be defined elsewhere by
            the subclass.

        regional : bool or None, optional
            Whether this is a regional or global grid.  If not provided here,
            it should be defined elsewhere by the subclass.
        """
        self.meshName = meshName
        self.regional = regional
        self.dims = None
        self.dimSizes = None
        self.coords = None

    def to_scrip(self, scripFileName):
        """
        Subclasses should overload this method to write a SCRIP file based on
        the mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written
        """
        raise NotImplemented('to_scrip is not implemented for this descriptor')

    def to_esmf(self, esmfFileName):
        """
        Subclasses should overload this method to write an ESMF mesh file based
        on the mesh.

        Parameters
        ----------
        esmfFileName : str
            The path to which the ESMF mesh file should be written
        """
        raise NotImplemented('to_esmf is not implemented for this descriptor')
