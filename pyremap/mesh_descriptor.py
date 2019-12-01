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
import xarray
import sys

from pyremap import write_netcdf


class MeshDescriptor(object):
    """
    A class for describing a mesh

    Attributes
    ----------
    meshname : str
        The name of the mesh, used as an attribute in NetCDF files and in
        mapping file names.

    lat_vertices, lon_vertices : numpy array
        latitude and longitude of vertices in radians, 1D arrays for meshes and
        2D for structured grids.

    x_vertices, y_vertices, z_vertices : numpy array
        Cartesian coordinates of vertices, 1D array for meshes and 2D for
        structured grids.

    vertices_on_cell : numpy array
        the connectivity describing which vertex indices surround each
        cell of a mesh.

    vertex_count_on_cells : numpy array
        The number of vertices on each cell for meshes.

    lat_cells, lon_cells : numpy array
        latitude and longitude of cells in radians, 1D arrays for meshes and
        2D for structured grids.

    x_cells, y_cells, z_cells : numpy array
        Cartesian coordinates of cells, 1D array for meshes and 2D for
        structured grids.

    cell_count : numpy array
        The number of cells, 1D for meshes and 2D for structured grids

    regional : bool
        If this is a regional or global grid

    sizes : OrderedDict
        A dictionary of cell dimensions and their sizes (1 entry for meshes,
        2 for grids)

    coords : OrderedDict
        A dictionary of xarray-style coordinates used to supply destination
        datasets and NetCDF files with appropriate coordinates after remapping

    history : str
        A history string to be written to output files for provenance.  When
        grids are read from files that contain a history attribute, the previous
        history should be prepended on this string.
    """

    def __init__(self):
        """
        Constructor simply sets all attributes to None.
        """

        self.meshname = None
        # variables related to vertices on cells of a grid or mesh
        self.lat_vertices = None
        self.lon_vertices = None
        self.x_vertices = None
        self.y_vertices = None
        self.z_vertices = None
        # variables related to vertices on cells of a mesh
        self.vertices_on_cell = None
        self.vertex_count_on_cells = None

        # variables related to cell centers of a grid or mesh
        self.lat_cells = None
        self.lon_cells = None
        self.x_cells = None
        self.y_cells = None
        self.z_cells = None

        self.cell_count = None

        self.regional = None

        self.sizes = None
        self.coords = None

        self.history = ' '.join(sys.argv[:])

    def set_lon_lat_vertices(self, lon_vertices, lat_vertices,
                             vertices_on_cell=None,
                             vertex_count_on_cells=None,
                             degrees=False):
        """
        Set the latitude and longitude of vertices around each finite volume
        cell.
        Parameters
        ----------
        lon_vertices, lat_vertices : numpy array
            latitude and longitude of vertices.  1D array for meshes, 2D for
            structured grids.  Each array should have the same size.  Use
            ``numpy.meshgrid`` or ``numpy.ndgrid`` to turn 1D grid dimensions
            into 2D arrays.

        vertices_on_cell : numpy array, optional
            the connectivity describing which vertex indices surround each
            cell of a mesh.  The first dimension is the number of cells and
            the second dimension is the maximum number of vertices on a cell.
            For grids, connectivity is determined automatically from the 2D
            array.

        vertex_count_on_cells : numpy array, optional
            The number of vertices on each cell for meshes.  If
            ``vertices_on_cell`` is provided but this variable is not, all
            cells are assumed to have the same number of vertices equal to the
            second dimension of ``vertices_on_cell``.

        degrees : bool, optional
            Whether latitude and longitude are given in degrees (as opposed to
            radians)
        """

        if lat_vertices.shape != lon_vertices.shape:
            raise ValueError('lon_vertices and lat_vertices must have the same '
                             'shape')

        self._set_cell_count_from_vertices(vertices_on_cell, lat_vertices)

        if degrees:
            self.lon_vertices = numpy.deg2rad(lon_vertices)
            self.lat_vertices = numpy.deg2rad(lat_vertices)
        else:
            self.lon_vertices = lon_vertices
            self.lat_vertices = lat_vertices

        self.vertices_on_cell = vertices_on_cell
        self.vertex_count_on_cells = vertex_count_on_cells

    def set_cartesian_vertices(self, x_vertices, y_vertices, z_vertices,
                               vertices_on_cell=None,
                               vertex_count_on_cells=None):
        """
        Set the Cartesian (x, y, z) coordinates of vertices around each finite
        volume cell.
        Parameters
        ----------
        x_vertices, y_vertices, z_vertices : numpy array
            Cartesian coordinates of vertices.  1D array for meshes, 2D for
            structured grids.  Each array should have the same size.  Use
            ``numpy.meshgrid`` or ``numpy.ndgrid`` to turn 1D grid dimensions
            into 2D arrays.

        vertices_on_cell : numpy array, optional
            the connectivity describing which vertex indices surround each
            cell of a mesh.  The first dimension is the number of cells and
            the second dimension is the maximum number of vertices on a cell.
            For grids, connectivity is determined automatically from the 2D
            array.

        vertex_count_on_cells : numpy array, optional
            The number of vertices on each cell for meshes.  If
            ``vertices_on_cell`` is provided but this variable is not, all
            cells are assumed to have the same number of vertices equal to the
            second dimension of ``vertices_on_cell``.
        """
        if x_vertices.shape != y_vertices.shape or \
                x_vertices.shape != z_vertices.shape:
            raise ValueError('x_vertices, y_vertices and z_vertices must have '
                             'the same shape')

        self._set_cell_count_from_vertices(vertices_on_cell, x_vertices)

        self.x_vertices, self.y_vertices, self.z_vertices = \
            MeshDescriptor._normalize(x_vertices, y_vertices, z_vertices)

        self.vertices_on_cell = vertices_on_cell
        self.vertex_count_on_cells = vertex_count_on_cells

    def set_lon_lat_cell_centers(self, lon_cells, lat_cells, degrees=False):
        """
        Set the latitude and longitude of vertices around each finite volume
        cell.
        Parameters
        ----------
        lon_cells, lat_cells : numpy array
            latitude and longitude of cell centers.  1D array for meshes, 2D for
            structured grids.  Each array should have the same size.  Use
            ``numpy.meshgrid`` or ``numpy.ndgrid`` to turn 1D grid dimensions
            into 2D arrays.

        degrees : bool, optional
            Whether latitude and longitude are given in degrees (as opposed to
            radians)
        """

        if lat_cells.shape != lon_cells.shape:
            raise ValueError('lat_cells and lon_cells must have the same '
                             'shape')

        self._set_cell_count_from_cells(lat_cells)

        if degrees:
            self.lon_cells = numpy.deg2rad(lon_cells)
            self.lat_cells = numpy.deg2rad(lat_cells)
        else:
            self.lon_cells = lon_cells
            self.lat_cells = lat_cells

    def set_cartesian_cell_centers(self, x_cells, y_cells, z_cells):
        """
        Set the Cartesian (x, y, z) coordinates of vertices around each finite
        volume cell.
        Parameters
        ----------
        x_cells, y_cells, z_cells : numpy array
            Cartesian coordinates of vertices.  1D array for meshes, 2D for
            structured grids.  Each array should have the same size.  Use
            ``numpy.meshgrid`` or ``numpy.ndgrid`` to turn 1D grid dimensions
            into 2D arrays.
        """
        if x_cells.shape != y_cells.shape or \
                x_cells.shape != z_cells.shape:
            raise ValueError('x_cells, y_cells and z_cells must have '
                             'the same shape')

        self._set_cell_count_from_cells(x_cells)

        self.x_cells, self.y_cells, self.z_cells = \
            MeshDescriptor._normalize(x_cells, y_cells, z_cells)

    def to_scrip(self, scrip_file_name):
        """
        Write out a SCRIP file that can be used to create a mapping file
        with ESMF_RegridWeightGen using this mesh as a source or destination.

        Parameters
        ----------
        scrip_file_name : str
            The path to which the SCRIP file should be written
        """
        if self.cell_count is None:
            raise ValueError('No coordinates have been set.')

        self._try_cartesian_to_lat_lon()

        ds = xarray.Dataset()

        ds['grid_dims'] = (('grid_rank',), self.cell_count[::-1])

        if self.lat_cells is not None and self.lon_cells is not None:
            ds['grid_center_lat'] = (('grid_size',), self.lat_cells.flat)
            ds.grid_center_lat.attrs['units'] = 'radians'
            ds['grid_center_lon'] = (('grid_size',), self.lon_cells.flat)
            ds.grid_center_lon.attrs['units'] = 'radians'

        if self.lat_vertices is not None and self.lon_vertices is not None:
            ds['grid_corner_lat'] = (('grid_size', 'grid_corners'),
                                     self._unravel_vertices(self.lat_vertices))
            ds.grid_corner_lat.attrs['units'] = 'radians'
            ds['grid_corner_lon'] = (('grid_size', 'grid_corners'),
                                     self._unravel_vertices(self.lon_vertices))
            ds.grid_corner_lon.attrs['units'] = 'radians'

        cell_count = numpy.product(self.cell_count)
        ds['grid_imask'] = (('grid_size',), numpy.ones(cell_count, dtype=int))
        ds.grid_imask.attrs['units'] = 'unitless'

        if self.meshname is not None:
            ds.attrs['meshname'] = self.meshname

        ds.attrs['history'] = self.history
        write_netcdf(ds, scrip_file_name, format='NETCDF4')

    def _set_cell_count_from_vertices(self, vertices_on_cell,
                                      field_on_vertices):
        if vertices_on_cell is None:
            # this is a grid, so lat_vertices should be 2D and the number of
            # cells in each dimension should be one less than the number of
            # vertices
            cell_count = numpy.array(field_on_vertices.shape) - 1
            if len(cell_count) != 2:
                raise ValueError('Expected 2D arrays for vertex coordinates for'
                                 ' a structured grid')
            if self.cell_count is None:
                self.cell_count = cell_count
            else:
                if numpy.any(cell_count != self.cell_count):
                    raise ValueError('Size of vertex coordinates should be 1 '
                                     'more in each dimension than the number '
                                     'of cells.')
        else:
            # this is a mesh, so lat_vertices should be 1D
            cell_count = numpy.array([vertices_on_cell.shape[0]])
            if self.cell_count is None:
                self.cell_count = cell_count
            else:
                if numpy.any(cell_count != self.cell_count):
                    raise ValueError('The first dimension of vertices_on_cell '
                                     'should be the number of cells.')

    def _set_cell_count_from_cells(self, field_on_cells):
        cell_count = numpy.array(field_on_cells.shape)
        if self.cell_count is None:
            self.cell_count = cell_count
        else:
            if numpy.any(cell_count != self.cell_count):
                raise ValueError('Size of cell coordinates do not match number'
                                 'of cells inferred from previously set '
                                 'coordinates.')

    def _try_cartesian_to_lat_lon(self):
        if self.x_vertices is not None and \
                self.y_vertices is not None and \
                self.z_vertices is not None:
            if self.lat_vertices is None or self.lon_vertices is None:
                self.lat_vertices, self.lon_vertices = \
                    MeshDescriptor._cartesian_to_lat_lon(
                        self.x_vertices, self.y_vertices, self.z_vertices)

        if self.x_cells is not None and \
                self.y_cells is not None and \
                self.z_cells is not None:
            if self.lat_cells is None or self.lon_cells is None:
                self.lat_cells, self.lon_cells = \
                    MeshDescriptor._cartesian_to_lat_lon(
                        self.x_cells, self.y_cells, self.z_cells)

    def _unravel_vertices(self, field_on_vertices):
        if self.vertices_on_cell is None:
            # this is a structured grid, so use _unwrap_corners to unravel
            unraveled = MeshDescriptor._unwrap_corners(field_on_vertices)
        elif self.vertex_count_on_cells is None:
            # this is an unstructured mesh and all cells have the same
            # number of vertices
            unraveled = field_on_vertices[self.vertices_on_cell]
        else:
            # this is an unstructured mesh, and polygons have different numbers
            # of vertices so we have to copy the last vertex to make them all
            # the same
            cell_count, max_verts = self.vertices_on_cell.shape
            unraveled = numpy.zeros(self.vertices_on_cell.shape)
            for vert_index in range(max_verts):
                cell_indices = numpy.arange(cell_count)
                # repeat last vertex where vert_index > vertex_count_on_cells
                local_vert_index = numpy.minimum(self.vertex_count_on_cells-1,
                                                 vert_index)
                vert_indices = \
                    self.vertices_on_cell[cell_indices, local_vert_index]
                unraveled[cell_indices, vert_index] = \
                    field_on_vertices[vert_indices]

        return unraveled

    @staticmethod
    def _cartesian_to_lat_lon(x, y, z):
        """Compute lat and lon in radians from cartesian coords."""
        lat = numpy.arcsin(z)
        lon = numpy.arctan2(y, x)
        return lat, lon

    @staticmethod
    def _normalize(x, y, z):
        """
        Normalize Cartesian coordinates to be on the unit sphere
        Parameters
        ----------
        x, y, z : numpy array
            The coordinates to normalize
        Returns
        -------
        x, y, z : numpy array
            The normalized coordinates on the unit sphere
        """
        magnitude = numpy.sqrt(x**2 + y**2 + z**2)
        # prevent division by zero
        magnitude[magnitude == 0.] = 1.
        x = x/magnitude
        y = y/magnitude
        z = z/magnitude
        return x, y, z

    @staticmethod
    def interp_extrap_corner(in_field):
        """
        Interpolate/extrapolate a 1D, regularly spaced field from grid cell
        centers to grid corners (vertices)
        Parameters
        ----------
        in_field : numpy array
            A field at cell centers in 1D

        Returns
        -------
        out_field : numpy array
            The field interpolated or extrapolated to 1D corners
        """
        out_field = numpy.zeros((len(in_field), 2))
        # interpolate the middle
        out_field[1:, 0] = 0.5 * (in_field[0:-1] + in_field[1:])
        out_field[0:-1, 1] = out_field[1:, 0]
        # extrapolate the ends
        out_field[0, 0] = 1.5 * in_field[0] - 0.5 * in_field[1]
        out_field[-1, 1] = 1.5 * in_field[-1] - 0.5 * in_field[-2]
        return out_field

    @staticmethod
    def interp_extrap_corners_2d(in_field):
        """
        Interpolate/extrapolate a 2D, regularly spaced field from grid cell
        centers to grid corners (vertices)
        Parameters
        ----------
        in_field : numpy array
            A field at cell centers in 2D

        Returns
        -------
        out_field : numpy array
            The field interpolated or extrapolated to corners
        """
        temp = numpy.zeros((in_field.shape[0], in_field.shape[1] + 1))
        temp[:, 1:-1] = 0.5 * (in_field[:, 0:-1] + in_field[:, 1:])
        # extrapolate the ends
        temp[:, 0] = 1.5 * in_field[:, 0] - 0.5 * in_field[:, 1]
        temp[:, -1] = 1.5 * in_field[:, -1] - 0.5 * in_field[:, -2]

        out_field = numpy.zeros((in_field.shape[0] + 1, in_field.shape[1] + 1))
        out_field[1:-1, :] = 0.5 * (temp[0:-1, :] + temp[1:, :])
        # extrapolate the ends
        out_field[0, :] = 1.5 * temp[0, :] - 0.5 * temp[1, :]
        out_field[-1, :] = 1.5 * temp[-1, :] - 0.5 * temp[-2, :]

        return out_field

    @staticmethod
    def _unwrap_corners(in_field):
        """
        Turn a 2D array of corners into an array of rectangular mesh elements
        """
        cell_count = (in_field.shape[0] - 1) * (in_field.shape[1] - 1)
        out_field = numpy.zeros((cell_count, 4))
        # corners are counterclockwise
        out_field[:, 0] = in_field[0:-1, 0:-1].flat
        out_field[:, 1] = in_field[0:-1, 1:].flat
        out_field[:, 2] = in_field[1:, 1:].flat
        out_field[:, 3] = in_field[1:, 0:-1].flat

        return out_field

    @staticmethod
    def _round_res(res):
        """Round the resolution to a reasonable number for grid names"""
        rounded = numpy.round(res*1000.)/1000.
        return '{}'.format(rounded)
