import os

from pyremap.remapper.build_map import _build_map
from pyremap.remapper.ncremap import _ncremap
from pyremap.remapper.remap_numpy import _remap_numpy
from pyremap.remapper.setup import _setup_remapper


class Remapper:
    """
    A class for remapping fields using a given mapping file.  The weights and
    indices from the mapping file can be loaded once and reused multiple times
    to map several fields between the same source and destination grids.

    Attributes
    ----------
    ntasks : int
        the target number of MPI tasks to use

    src_grid_info : dict
        Information about the source grid

    dst_grid_info : dict
        Information about the destination grid

    method : {'bilinear', 'neareststod', 'conserve'}
        The method of interpolation used

    map_filename : str or None
        The name of the output mapping file

    use_tmp : bool
        If True, use a temporary directory for the SCRIP files.  If False, the
        SCRIP files will be created in the current working directory.

    expand_dist : float or None
        The distance in meters over which to expand destination grid cells for
        smoothing

    expand_factor : float or None
        The factor by which to expand destination grid cells for smoothing

    src_scrip_filename : str
        The name of the SCRIP file for the source mesh

    dst_scrip_filename : str
        The name of the SCRIP file for the destination mesh

    format : {'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_DATA'}
        The NetCDF file format to use for scrip files.  The default is
        ``NETCDF3_64BIT_DATA`` for compatibility with MOAB.

    src_descriptor: pyremap.descriptor.MeshDescriptor
        The source mesh descriptor

    dst_descriptor: pyremap.descriptor.MeshDescriptor
        The destination mesh descriptor

    map_tool : {'esmf', 'moab'}
        The tool to use for building the mapping file.  The default is
        ``esmf``.

    esmf_path : str or None
        The path to the ESMF installation.  If None, ``ESMF_RegridWeigthGen``
        must be in the PATH.

    moab_path : str or None
        The path to the MOAB installation.  If None, ``mptempest`` and other
        MOAB tools must be in the PATH.

    parallel_exec : {'mpirun', 'srun'}
        The command to use for running the mapping tool.  The default is
        ``mpirun``.
    """  # noqa: E501

    def __init__(
        self,
        ntasks=1,
        map_filename=None,
        method='bilinear',
        use_tmp=True,
    ):
        """
        Create a new mapping object

        Parameters
        ----------
        ntasks : int
            the target number of MPI tasks to use

        map_filename : str, optional
            The name of the output mapping file,
            ``map_{source_type}_{dest_type}_{method}.nc`` by default

        method : {'bilinear', 'neareststod', 'conserve'}, optional
            The method of interpolation used

        use_tmp : bool, optional
            If True, use a temporary directory for the SCRIP files.  The
            default is False, which means the SCRIP files will be created in
            the current working directory.
        """  # noqa: E501
        self.ntasks = ntasks
        self.src_grid_info = dict()
        self.dst_grid_info = dict()
        self.map_filename = map_filename
        self.method = method
        self.use_tmp = use_tmp
        self.expand_dist = None
        self.expand_factor = None
        self.src_scrip_filename = 'src_mesh.nc'
        self.dst_scrip_filename = 'dst_mesh.nc'
        self.format = 'NETCDF3_64BIT_DATA'
        self.src_descriptor = None
        self.dst_descriptor = None
        self.map_tool = 'esmf'
        self.esmf_path = None
        self.moab_path = None
        self.parallel_exec = 'mpirun'
        self._ds_map = None
        self._matrix = None

    def src_from_lon_lat(
        self, filename, mesh_name=None, lon_var='lon', lat_var='lat'
    ):
        """
        Set the source grid from a file with a longitude-latitude grid.  The
        latitude and longitude variables can be 1D or 2D.

        Parameters
        ----------
        filename : str
            A file containing the latitude-longitude grid

        mesh_name : str, optional
            The name of the lon-lat grid (defaults to resolution and units,
            something like "0.5x0.5degree")

        lon_var : str, optional
            The name of the longitude coordinate in the file

        lat_var : str, optional
            The name of the latitude coordinate in the file
        """
        src = dict()
        src['type'] = 'lon-lat'
        src['filename'] = filename
        src['lon'] = lon_var
        src['lat'] = lat_var
        if mesh_name is not None:
            src['name'] = mesh_name
        self.src_grid_info = src

    def dst_from_lon_lat(
        self, filename, mesh_name=None, lon_var='lon', lat_var='lat'
    ):
        """
        Set the destination grid from a file with a longitude-latitude grid.
        The latitude and longitude variables can be 1D or 2D.

        Parameters
        ----------
        filename : str
            A file containing the latitude-longitude grid

        mesh_name : str, optional
            The name of the lon-lat grid (defaults to resolution and units,
            something like "0.5x0.5degree")

        lon_var : str, optional
            The name of the longitude coordinate in the file

        lat_var : str, optional
            The name of the latitude coordinate in the file
        """
        dst = dict()
        dst['type'] = 'lon-lat'
        dst['filename'] = filename
        dst['lon'] = lon_var
        dst['lat'] = lat_var
        if mesh_name is not None:
            dst['name'] = mesh_name
        self.dst_grid_info = dst

    def dst_global_lon_lat(self, dlon, dlat, lon_min=-180.0, mesh_name=None):
        """
        Set the destination grid from a file with a longitude-latitude grid.
        The latitude and longitude variables can be 1D or 2D.

        Parameters
        ----------
        dlon : float
            The longitude resolution in degrees

        dlat : float
            The latitude resolution in degrees

        lon_min : float, optional
            The longitude for the left-hand edge of the global grid in degrees

        mesh_name : str, optional
            The name of the lon-lat grid (defaults to resolution and units,
            something like "0.5x0.5degree")
        """

        dst = dict()
        dst['type'] = 'lon-lat'
        dst['dlon'] = dlon
        dst['dlat'] = dlat
        dst['lon_min'] = lon_min
        if mesh_name is not None:
            dst['name'] = mesh_name
        self.dst_grid_info = dst

    def src_from_proj(
        self,
        filename,
        mesh_name,
        x_var='x',
        y_var='y',
        proj_attr=None,
        proj_str=None,
    ):
        """
        Set the source grid from a file with a projection grid.

        Parameters
        ----------
        filename : str
            A file containing the projection grid

        mesh_name : str
            The name of the projection grid

        x_var : str, optional
            The name of the x coordinate in the file

        y_var : str, optional
            The name of the y coordinate in the file

        proj_attr : str, optional
            The name of a global attribute in the file containing the proj
            string for the projection

        proj_str : str, optional
            A proj string defining the projection, ignored if ``proj_attr``
            is provided
        """
        src = dict()
        src['type'] = 'proj'
        src['filename'] = filename
        src['name'] = mesh_name
        src['x'] = x_var
        src['y'] = y_var
        if proj_attr is not None:
            src['proj_attr'] = proj_attr
        elif proj_str is not None:
            src['proj_str'] = proj_str
        else:
            raise ValueError('Must provide one of "proj_attr" or "proj_str".')
        self.src_grid_info = src

    def dst_from_proj(
        self,
        filename,
        mesh_name,
        x_var='x',
        y_var='y',
        proj_attr=None,
        proj_str=None,
    ):
        """
        Set the destination grid from a file with a projection grid.

        Parameters
        ----------
        filename : str
            A file containing the projection grid

        mesh_name : str
            The name of the projection grid

        x_var : str, optional
            The name of the x coordinate in the file

        y_var : str, optional
            The name of the y coordinate in the file

        proj_attr : str, optional
            The name of a global attribute in the file containing the proj
            string for the projection

        proj_str : str, optional
            A proj string defining the projection, ignored if ``proj_attr``
            is provided
        """
        dst = dict()
        dst['type'] = 'proj'
        dst['filename'] = filename
        dst['name'] = mesh_name
        dst['x'] = x_var
        dst['y'] = y_var
        if proj_attr is not None:
            dst['proj_attr'] = proj_attr
        elif proj_str is not None:
            dst['proj_str'] = proj_str
        else:
            raise ValueError('Must provide one of "proj_attr" or "proj_str".')
        self.dst_grid_info = dst

    def dst_from_points(
        self, filename, mesh_name, lon_var='lon', lat_var='lat'
    ):
        """
        Set the destination grid from a file with a collection of points.

        Parameters
        ----------
        filename : str
            A file containing the latitude-longitude grid

        mesh_name : str
            The name of the point collection

        lon_var : str, optional
            The name of the longitude coordinate in the file

        lat_var : str, optional
            The name of the latitude coordinate in the file
        """
        dst = dict()
        dst['type'] = 'points'
        dst['filename'] = filename
        dst['name'] = mesh_name
        dst['lon'] = lon_var
        dst['lat'] = lat_var
        self.dst_grid_info = dst

    def src_from_mpas(self, filename, mesh_name, mesh_type='cell'):
        """
        Set the source grid from an MPAS mesh file

        Parameters
        ----------
        filename : str
            A file containing the MPAS mesh

        mesh_name : str
            The name of the MPAS mesh

        mesh_type : {'cell', 'edge', 'vertex'}, optional
            Which type of MPAS mesh
        """
        src = dict()
        src['type'] = 'mpas'
        src['filename'] = filename
        src['name'] = mesh_name
        src['mpas_mesh_type'] = mesh_type
        self.src_grid_info = src

    def dst_from_mpas(self, filename, mesh_name, mesh_type='cell'):
        """
        Set the destination grid from an MPAS mesh file

        Parameters
        ----------
        filename : str
            A file containing the MPAS mesh

        mesh_name : str
            The name of the MPAS mesh

        mesh_type : {'cell', 'edge', 'vertex'}, optional
            Which type of MPAS mesh
        """
        dst = dict()
        dst['type'] = 'mpas'
        dst['filename'] = filename
        dst['name'] = mesh_name
        dst['mpas_mesh_type'] = mesh_type
        self.dst_grid_info = dst

    def build_map(self):
        """
        Make the mapping file
        """
        _build_map(self)

    def ncremap(
            self,
            in_filename,
            out_filename,
            variable_list=None,
            overwrite=False,
            renormalize=None,
            logger=None,
            replace_mpas_fill=False,
            parallel_exec=None):
        """
        Given a source file defining either an MPAS mesh or a lat-lon grid and
        a destination file or set of arrays defining a lat-lon grid, constructs
        a mapping file used for interpolation between the source and
        destination grids.

        Parameters
        ----------
        in_filename : str
            The path to the file containing a data set on the source grid

        out_filename : str
            The path where the data on the destination grid should be written

        variable_list : list of str, optional
            A list of variables to be mapped.  By default, all variables are
            mapped

        overwrite : bool, optional
            Whether the destination file should be overwritten if it already
            exists. If `False`, and the destination file is already present,
            the function does nothing and returns immediately

        renormalize : float, optional
            A threshold to use to renormalize the data

        logger : ``logging.Logger``, optional
            A logger to which ncclimo output should be redirected

        replace_mpas_fill : bool, optional
            For MPAS meshes, whether add a ``_FillValue`` attribute (missing
            from MPAS output).  If this has been handled before the call,
            replacing the fill value again may cause errors.

        parallel_exec : {'srun'}, optional
            The name of the parallel executable to use to launch ncremap.
            By default, none is used.

        Raises
        ------
        OSError
            If ``ncremap`` is not in the system path.

        ValueError
            If ``mapping_file_name`` is ``None`` (meaning no remapping is
            needed).
        """
        _setup_remapper(self)
        if not os.path.exists(self.map_filename):
            _build_map(self)

        _ncremap(
            self,
            in_filename,
            out_filename,
            variable_list,
            overwrite,
            renormalize,
            logger,
            replace_mpas_fill,
            parallel_exec
        )

    def remap_numpy(self, ds, renormalization_threshold=None):
        """
        Given a source data set, returns a remapped version of the data set,
        possibly masked and renormalized.

        Parameters
        ----------
        ds : xarray.Dataset or xarray.DataArray
            The dimention(s) along ``self.sourceDimNames`` must match
            ``self.src_grid_dims`` read from the mapping file.

        renormalization_threshold : float, optional
            The minimum weight of a denstination cell after remapping, below
            which it is masked out, or ``None`` for no renormalization and
            masking.

        Returns
        -------
        ds_remap : xarray.Dataset or xarray.DataArray
            Returns a remapped data set (or data array) where dimensions other
            than ``self.src_descriptor.dims`` are the same as in ``ds`` and the
            dimension(s) given by ``self.src_descriptor.dims`` have been
            replaced by ``self.dst_descriptor.dims``.
        """
        return _remap_numpy(self, ds, renormalization_threshold)
