# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2025 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2025 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2025 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/pyremap/main/LICENSE

from typing import Optional

import numpy as np
import pyproj
import xarray as xr

from pyremap.descriptor.mesh_descriptor import MeshDescriptor
from pyremap.descriptor.utility import (
    add_history,
    expand_scrip,
    interp_extrap_corner,
    unwrap_corners,
)


class ProjectionGridDescriptor(MeshDescriptor):
    """
    A class for describing a general logically rectangular grid that can be
    defined by a ``pyproj`` projection.

    Attributes
    ----------
    projection : pyproj.Proj
        The projection used to map from grid x-y space to latitude and
        longitude

    lat_lon_projection : pyproj.Proj
        lat-lon projection used to transform from x-y to lat-lon space

    x : numpy.ndarray
        The latitude coordinate at grid-cell centers

    y : numpy.ndarray
        The longitude coordinate at grid-cell centers

    x_corner : numpy.ndarray
        The latitude coordinate at grid-cell corners

    y_corner : numpy.ndarray
        The longitude coordinate at grid-cell corners

    history : str
        The history attribute written to SCRIP files

    x_var_name : str
        The name of the x variable

    y_var_name : str
        The name of the y variable
    """

    def __init__(self, projection, mesh_name=None):
        """
        Constructor stores the projection

        Parameters
        ----------
        projection : pyproj.Proj
            The projection used to map from grid x-y space to latitude and
            longitude

        mesh_name : str, optional
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``)
        """
        super().__init__(mesh_name=mesh_name, regional=True)

        self.projection = projection
        self.lat_lon_projection = pyproj.Proj(proj='latlong', datum='WGS84')

        self.x: Optional[np.ndarray] = None
        self.y: Optional[np.ndarray] = None
        self.x_corner: Optional[np.ndarray] = None
        self.y_corner: Optional[np.ndarray] = None
        self.history: Optional[str] = None
        self.x_var_name: Optional[str] = None
        self.y_var_name: Optional[str] = None

    @classmethod
    def read(
        cls,
        projection,
        filename,
        mesh_name=None,
        x_var_name='x',
        y_var_name='y',
    ):
        """
        Given a grid file with x and y coordinates defining the axes of the
        logically rectangular grid, read in the x and y coordinates and
        interpolate/extrapolate to locate corners.

        Parameters
        ----------
        projection : pyproj.Proj object
            The projection used to map from grid x-y space to latitude and
            longitude

        filename : str
            The path of the file containing the grid data

        mesh_name : str, optional
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``).  If not
            provided, the data set in ``filename`` must have a global
            attribute ``mesh_name`` of ``meshName`` that will be used instead.

        x_var_name, y_var_name : str, optional
            The name of the x and y (in meters) variables in the grid file
        """
        ds = xr.open_dataset(filename)

        descriptor = cls(projection, mesh_name=mesh_name)
        descriptor.mesh_name_from_attr(ds)
        if descriptor.mesh_name is None:
            raise ValueError('No mesh_name provided or found in file.')

        # Get info from input file
        descriptor.x = np.array(ds[x_var_name].values, float)
        descriptor.y = np.array(ds[y_var_name].values, float)

        descriptor._set_coords(
            x_var_name,
            y_var_name,
            ds[x_var_name].dims[0],
            ds[y_var_name].dims[0],
        )

        # interp/extrap corners
        descriptor.x_corner = interp_extrap_corner(descriptor.x)
        descriptor.y_corner = interp_extrap_corner(descriptor.y)

        descriptor.history = add_history(ds=ds)
        return descriptor

    @classmethod
    def create(cls, projection, x, y, mesh_name):
        """
        Given x and y coordinates defining the axes of the logically
        rectangular grid, save the coordinates interpolate/extrapolate to
        locate corners.

        Parameters
        ----------
        projection : pyproj.Proj object
            The projection used to map from grid x-y space to latitude and
            longitude

        x : numpy.ndarray
            One dimensional array defining the x coordinate of grid cell
            centers.

        y : numpy.ndarray
            One dimensional array defining the y coordinate of grid cell
            centers.

        mesh_name : str
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``)
        """
        descriptor = cls(projection, mesh_name=mesh_name)

        descriptor.x = x
        descriptor.y = y

        descriptor._set_coords('x', 'y', 'x', 'y')

        # interp/extrap corners
        descriptor.x_corner = interp_extrap_corner(descriptor.x)
        descriptor.y_corner = interp_extrap_corner(descriptor.y)
        descriptor.history = add_history()
        return descriptor

    def to_scrip(self, scrip_filename, expand_dist=None, expand_factor=None):
        """
        Create a SCRIP file based on the grid and projection.

        Parameters
        ----------
        scrip_filename : str
            The path to which the SCRIP file should be written

        expand_dist : float or numpy.ndarray, optional
            A distance in meters to expand each grid cell outward from the
            center.  If a ``numpy.ndarray``, one value per cell.

        expand_factor : float or numpy.ndarray, optional
            A factor by which to expand each grid cell outward from the center.
            If a ``numpy.ndarray``, one value per cell.
        """
        assert self.x is not None, 'x must be set before calling to_scrip'
        assert self.y is not None, 'y must be set before calling to_scrip'
        assert self.x_corner is not None, (
            'x_corner must be set before calling to_scrip'
        )
        assert self.y_corner is not None, (
            'y_corner must be set before calling to_scrip'
        )

        ds = xr.Dataset()

        nx = len(self.x)
        ny = len(self.y)

        grid_size = nx * ny

        (center_x, center_y) = np.meshgrid(self.x, self.y)
        (corner_x, corner_y) = np.meshgrid(self.x_corner, self.y_corner)

        (center_lat, center_lon) = self.project_to_lat_lon(center_x, center_y)
        (corner_lat, corner_lon) = self.project_to_lat_lon(corner_x, corner_y)

        ds['grid_center_lat'] = (('grid_size',), center_lat.flat)
        ds['grid_center_lon'] = (('grid_size',), center_lon.flat)
        ds['grid_corner_lat'] = (
            ('grid_size', 'grid_corners'),
            unwrap_corners(corner_lat),
        )
        ds['grid_corner_lon'] = (
            ('grid_size', 'grid_corners'),
            unwrap_corners(corner_lon),
        )

        ds['grid_dims'] = xr.DataArray([nx, ny], dims=('grid_rank',)).astype(
            'int32'
        )

        ds['grid_imask'] = xr.DataArray(
            np.ones(grid_size, dtype='int32'), dims=('grid_size',)
        )

        if expand_dist is not None or expand_factor is not None:
            expand_scrip(ds, expand_dist, expand_factor)

        ds.grid_center_lat.attrs['units'] = 'degrees'
        ds.grid_center_lon.attrs['units'] = 'degrees'
        ds.grid_corner_lat.attrs['units'] = 'degrees'
        ds.grid_corner_lon.attrs['units'] = 'degrees'
        ds.grid_imask.attrs['units'] = 'unitless'

        ds.attrs['mesh_name'] = self.mesh_name
        ds.attrs['history'] = self.history
        self.write_netcdf(ds, scrip_filename)

    def project_to_lat_lon(self, x, y):
        """
        Given x and y locations of points in a projection, returns the
        corresponding latitude and longitude of each point.

        Parameters
        ----------
        x : numpy.ndarray
            x array of points in the projection

        y : numpy.ndarray
            y array of points in the projection

        Returns
        -------
        lat : numpy.ndarray
            The latitude in degrees with the same size as x and y

        lon : numpy.ndarray
            The longitude in degrees with the same size as x and y
        """
        transformer = pyproj.Transformer.from_proj(
            self.projection, self.lat_lon_projection
        )
        lon, lat = transformer.transform(x, y)

        return lat, lon

    def _set_coords(self, x_var_name, y_var_name, x_dim_name, y_dim_name):
        """
        Set up a coords dict with x, y, lat and lon
        """
        assert self.x is not None, 'x must be set before calling to_scrip'
        assert self.y is not None, 'y must be set before calling to_scrip'
        self.x_var_name = x_var_name
        self.y_var_name = y_var_name
        (x, y) = np.meshgrid(self.x, self.y)
        (lat, lon) = self.project_to_lat_lon(x, y)

        self.coords = {
            x_var_name: {
                'dims': x_dim_name,
                'data': self.x,
                'attrs': {'units': 'meters'},
            },
            y_var_name: {
                'dims': y_dim_name,
                'data': self.y,
                'attrs': {'units': 'meters'},
            },
            'lat': {
                'dims': (y_dim_name, x_dim_name),
                'data': lat,
                'attrs': {'units': 'degrees'},
            },
            'lon': {
                'dims': (y_dim_name, x_dim_name),
                'data': lon,
                'attrs': {'units': 'degrees'},
            },
        }

        self.dims = [y_dim_name, x_dim_name]
        self.dim_sizes = [len(self.y), len(self.x)]
