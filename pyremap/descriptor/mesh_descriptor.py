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

from pyremap.utility import write_netcdf as write_netcdf_util


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
              'NETCDF3_CLASSIC'}
        The NetCDF file format to use.  Default is ``'NETCDF4'``

    engine : {'netcdf4', 'scipy', 'h5netcdf'}
        The library to use for xarray NetCDF output.  The default is
        ``'netcdf4'``

    logger : logging.Logger or None
        A logger for command-line output.  If None, the logger will be set to
        stdout/stderr.
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
        self.logger = None

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

    def write_netcdf(self, ds, filename):
        """
        Write the mesh to a NetCDF file

        Parameters
        ----------
        ds : xarray.Dataset
            The dataset to save

        filename : str
            The path for the NetCDF file to write

        """
        write_netcdf_util(
            ds,
            filename,
            format=self.format,
            engine=self.engine,
            logger=self.logger
        )
