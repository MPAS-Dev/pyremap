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

    format : {'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT',
              'NETCDF3_64BIT_OFFSET', 'NETCDF3_64BIT_DATA', 'NETCDF3_CLASSIC'}
        The NetCDF file format to use.  Default is ``'NETCDF4'``

    engine : {'netcdf4', 'scipy', 'h5netcdf'}
        The library to use for xarray NetCDF output.  The default is
        ``'netcdf4'``
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
        self.format = 'NETCDF4'
        self.engine = None

    def to_scrip(self, scripFileName, expandDist=None, expandFactor=None):
        """
        Subclasses should overload this method to write a SCRIP file based on
        the mesh.

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
        raise NotImplementedError(
            'to_scrip is not implemented for this descriptor')

    def to_esmf(self, esmfFileName):
        """
        Subclasses should overload this method to write an ESMF mesh file based
        on the mesh.

        Parameters
        ----------
        esmfFileName : str
            The path to which the ESMF mesh file should be written
        """
        raise NotImplementedError(
            'to_esmf is not implemented for this descriptor')
