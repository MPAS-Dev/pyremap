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

import sys

import netCDF4
import numpy as np
import pyproj.enums
from pyproj import Transformer


def interp_extrap_corner(inField):
    """Interpolate/extrapolate a 1D field from grid centers to grid corners"""

    outField = np.zeros(len(inField) + 1)
    outField[1:-1] = 0.5 * (inField[0:-1] + inField[1:])
    # extrapolate the ends
    outField[0] = 1.5 * inField[0] - 0.5 * inField[1]
    outField[-1] = 1.5 * inField[-1] - 0.5 * inField[-2]
    return outField


def interp_extrap_corners_2d(inField):
    """Interpolate/extrapolate a 1D field from grid centers to grid corners"""

    temp = np.zeros((inField.shape[0], inField.shape[1] + 1))
    temp[:, 1:-1] = 0.5 * (inField[:, 0:-1] + inField[:, 1:])
    # extrapolate the ends
    temp[:, 0] = 1.5 * inField[:, 0] - 0.5 * inField[:, 1]
    temp[:, -1] = 1.5 * inField[:, -1] - 0.5 * inField[:, -2]

    outField = np.zeros((inField.shape[0] + 1, inField.shape[1] + 1))
    outField[1:-1, :] = 0.5 * (temp[0:-1, :] + temp[1:, :])
    # extrapolate the ends
    outField[0, :] = 1.5 * temp[0, :] - 0.5 * temp[1, :]
    outField[-1, :] = 1.5 * temp[-1, :] - 0.5 * temp[-2, :]

    return outField


def expand_scrip(ds, expandDist, expandFactor):
    """
    Expand the vertices of cells outward from the center of the cell

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset for a SCRIP file

    expandDist : float or np.ndarray
        A distance in meters to expand each grid cell outward from the
        center.  If a ``np.ndarray``, one value per cell.

    expandFactor : float or np.ndarray
        A factor by which to expand each grid cell outward from the center.
        If a ``np.ndarray``, one value per cell.
    """
    grid_center_lat = ds.grid_center_lat
    grid_center_lon = ds.grid_center_lon
    grid_corner_lat = ds.grid_corner_lat
    grid_corner_lon = ds.grid_corner_lon
    grid_size = ds.sizes['grid_size']
    grid_corners = ds.sizes['grid_corners']

    radians = 'rad' in grid_center_lat.attrs['units']

    trans_lon_lat_to_xyz = Transformer.from_crs(4979, 4978, always_xy=True)
    x_center, y_center, z_center = trans_lon_lat_to_xyz.transform(
        grid_center_lon.values,
        grid_center_lat.values,
        np.zeros(grid_size),
        radians=radians)

    x_corner, y_corner, z_corner = trans_lon_lat_to_xyz.transform(
        grid_corner_lon.values,
        grid_corner_lat.values,
        np.zeros((grid_size, grid_corners)),
        radians=radians)

    if expandFactor is None:
        expandFactor = 1.

    if expandDist is None:
        expandDist = 0.

    for index in range(grid_corners):
        dx = x_corner[:, index] - x_center
        dy = y_corner[:, index] - y_center
        dz = z_corner[:, index] - z_center
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        factor = (expandFactor * dist + expandDist) / dist
        x_corner[:, index] = factor * dx + x_center
        y_corner[:, index] = factor * dy + y_center
        z_corner[:, index] = factor * dz + z_center

    grid_corner_lon, grid_corner_lat, _ = trans_lon_lat_to_xyz.transform(
        x_corner,
        y_corner,
        z_corner,
        radians=radians,
        direction=pyproj.enums.TransformDirection.INVERSE)

    ds['grid_corner_lat'] = (('grid_size', 'grid_corners'), grid_corner_lat)
    ds['grid_corner_lon'] = (('grid_size', 'grid_corners'), grid_corner_lon)


def unwrap_corners(inField):
    """Turn a 2D array of corners into an array of rectangular mesh elements"""
    outField = np.zeros(((inField.shape[0] - 1) * (inField.shape[1] - 1), 4))
    # corners are counterclockwise
    outField[:, 0] = inField[0:-1, 0:-1].flat
    outField[:, 1] = inField[0:-1, 1:].flat
    outField[:, 2] = inField[1:, 1:].flat
    outField[:, 3] = inField[1:, 0:-1].flat

    return outField


def round_res(res):
    """Round the resolution to a reasonable number for grid names"""
    rounded = np.round(res * 1000.) / 1000.
    return '{}'.format(rounded)


def add_history(ds=None):
    """Get the history attribute, possibly adding it to existing history"""
    history = ' '.join(sys.argv[:])
    if ds is not None and 'history' in ds.attrs:
        prev_hist = ds.attrs['history']
        if isinstance(prev_hist, np.ndarray):
            prev_hist = '\n'.join(prev_hist)
        history = '\n'.join([prev_hist, history])
    return history


def write_netcdf(ds, filename, format, engine=None, fillvalues=None):
    """
    Write an xarray.Dataset to a file with NetCDF4 fill values.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to save

    filename : str
        The path for the NetCDF file to write

    format : {'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC'}
        The NetCDF file format to use.  Default is
        ``mpas_tools.io.default_format``, which can be modified but which
        defaults to ``'NETCDF3_64BIT'``

    engine : {'netcdf4', 'scipy', 'h5netcdf'}, optional
        The library to use for NetCDF output.  The default is the same as
        in :py:meth:`xarray.Dataset.to_netcdf` and depends on ``format``.
        You can override the default by setting
        ``mpas_tools.io.default_engine``

    fillvalues : dict, optional
        A dictionary of fill values for different NetCDF types.  Default is
        ``mpas_tools.io.default_fills``, which can be modified but which
        defaults to ``netCDF4.default_fillvals``
    """

    if fillvalues is None:
        fillvalues = netCDF4.default_fillvals

    encoding_dict = {}
    variable_names = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variable_name in variable_names:
        is_numeric = np.issubdtype(ds[variable_name].dtype, np.number)
        if is_numeric and np.any(np.isnan(ds[variable_name])):
            dtype = ds[variable_name].dtype
            for fill_type in fillvalues:
                if dtype == np.dtype(fill_type):
                    encoding_dict[variable_name] = \
                        {'_FillValue': fillvalues[fill_type]}
                    break
        else:
            encoding_dict[variable_name] = {'_FillValue': None}

    ds.to_netcdf(
        filename, encoding=encoding_dict, format=format, engine=engine)
