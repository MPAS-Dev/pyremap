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

import sys

import numpy as np
import pyproj.enums
from pyproj import Transformer


def interp_extrap_corner(in_field):
    """Interpolate/extrapolate a 1D field from grid centers to grid corners"""

    out_field = np.zeros(len(in_field) + 1)
    out_field[1:-1] = 0.5 * (in_field[0:-1] + in_field[1:])
    # extrapolate the ends
    out_field[0] = 1.5 * in_field[0] - 0.5 * in_field[1]
    out_field[-1] = 1.5 * in_field[-1] - 0.5 * in_field[-2]
    return out_field


def interp_extrap_corners_2d(in_field):
    """Interpolate/extrapolate a 1D field from grid centers to grid corners"""

    temp = np.zeros((in_field.shape[0], in_field.shape[1] + 1))
    temp[:, 1:-1] = 0.5 * (in_field[:, 0:-1] + in_field[:, 1:])
    # extrapolate the ends
    temp[:, 0] = 1.5 * in_field[:, 0] - 0.5 * in_field[:, 1]
    temp[:, -1] = 1.5 * in_field[:, -1] - 0.5 * in_field[:, -2]

    out_field = np.zeros((in_field.shape[0] + 1, in_field.shape[1] + 1))
    out_field[1:-1, :] = 0.5 * (temp[0:-1, :] + temp[1:, :])
    # extrapolate the ends
    out_field[0, :] = 1.5 * temp[0, :] - 0.5 * temp[1, :]
    out_field[-1, :] = 1.5 * temp[-1, :] - 0.5 * temp[-2, :]

    return out_field


def expand_scrip(ds, expand_dist, expand_factor):
    """
    Expand the vertices of cells outward from the center of the cell

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset for a SCRIP file

    expand_dist : float or np.ndarray
        A distance in meters to expand each grid cell outward from the
        center.  If a ``np.ndarray``, one value per cell.

    expand_factor : float or np.ndarray
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
        radians=radians,
    )

    x_corner, y_corner, z_corner = trans_lon_lat_to_xyz.transform(
        grid_corner_lon.values,
        grid_corner_lat.values,
        np.zeros((grid_size, grid_corners)),
        radians=radians,
    )

    if expand_factor is None:
        expand_factor = 1.0

    if expand_dist is None:
        expand_dist = 0.0

    for index in range(grid_corners):
        dx = x_corner[:, index] - x_center
        dy = y_corner[:, index] - y_center
        dz = z_corner[:, index] - z_center
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        factor = (expand_factor * dist + expand_dist) / dist
        x_corner[:, index] = factor * dx + x_center
        y_corner[:, index] = factor * dy + y_center
        z_corner[:, index] = factor * dz + z_center

    grid_corner_lon, grid_corner_lat, _ = trans_lon_lat_to_xyz.transform(
        x_corner,
        y_corner,
        z_corner,
        radians=radians,
        direction=pyproj.enums.TransformDirection.INVERSE,
    )

    ds['grid_corner_lat'] = (('grid_size', 'grid_corners'), grid_corner_lat)
    ds['grid_corner_lon'] = (('grid_size', 'grid_corners'), grid_corner_lon)


def unwrap_corners(in_field):
    """Turn a 2D array of corners into an array of rectangular mesh elements"""
    out_field = np.zeros(
        ((in_field.shape[0] - 1) * (in_field.shape[1] - 1), 4)
    )
    # corners are counterclockwise
    out_field[:, 0] = in_field[0:-1, 0:-1].flat
    out_field[:, 1] = in_field[0:-1, 1:].flat
    out_field[:, 2] = in_field[1:, 1:].flat
    out_field[:, 3] = in_field[1:, 0:-1].flat

    return out_field


def round_res(res):
    """Round the resolution to a reasonable number for grid names"""
    rounded = np.round(res * 1000.0) / 1000.0
    return f'{rounded}'


def add_history(ds=None):
    """Get the history attribute, possibly adding it to existing history"""
    history = ' '.join(sys.argv[:])
    if ds is not None and 'history' in ds.attrs:
        prev_history = ds.attrs['history']
        if isinstance(prev_history, np.ndarray):
            prev_history = '\n'.join(prev_history)
        history = '\n'.join([prev_history, history])
    return history
