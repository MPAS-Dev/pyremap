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
import warnings

import numpy as np
import pyproj.enums
from pyproj import Transformer


def get_corners_1d(ds, var_name):
    """
    Get the coordinates of grid-cell corners along a 1D coordinate, using the
    CF ``bounds`` of the coordinate if they are available and usable, and
    interpolating and extrapolating from cell centers if not.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset containing the coordinate variable ``var_name``

    var_name : str
        The name of the 1D coordinate variable

    Returns
    -------
    corner : numpy.ndarray
        A 1D array of corner coordinates with one more element than
        ``ds[var_name]``
    """
    center = np.array(ds[var_name].values, float)
    bounds = _get_cf_bounds(ds, var_name, shape=(len(center), 2))
    if bounds is not None:
        corner = _corners_from_bounds_1d(bounds)
        if corner is not None:
            return corner
        warnings.warn(
            f'The CF bounds of {var_name} are not contiguous so corners '
            f'will be interpolated and extrapolated from cell centers '
            f'instead.',
            stacklevel=2,
        )

    return interp_extrap_corner(center)


def get_corners_2d(ds, lat_var_name, lon_var_name):
    """
    Get the coordinates of grid-cell corners for 2D latitude and longitude
    coordinates, using the CF ``bounds`` of the coordinates if they are
    available and usable, and interpolating and extrapolating from cell
    centers if not.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset containing the coordinate variables

    lat_var_name, lon_var_name : str
        The names of the 2D latitude and longitude variables

    Returns
    -------
    lat_corner, lon_corner : numpy.ndarray
        2D arrays of corner coordinates with one more element than
        ``ds[lat_var_name]`` and ``ds[lon_var_name]`` along each dimension
    """
    lat = np.array(ds[lat_var_name].values, float)
    lon = np.array(ds[lon_var_name].values, float)
    shape = (lat.shape[0], lat.shape[1], 4)
    lat_bounds = _get_cf_bounds(ds, lat_var_name, shape=shape)
    lon_bounds = _get_cf_bounds(ds, lon_var_name, shape=shape)
    if lat_bounds is not None and lon_bounds is not None:
        corners = _corners_from_bounds_2d(lat_bounds, lon_bounds)
        if corners is not None:
            return corners
        warnings.warn(
            f'The CF bounds of {lat_var_name} and {lon_var_name} do not '
            f'share vertices between neighboring cells so corners will be '
            f'interpolated and extrapolated from cell centers instead.',
            stacklevel=2,
        )
    elif lat_bounds is not None or lon_bounds is not None:
        warnings.warn(
            f'Only one of {lat_var_name} and {lon_var_name} has usable CF '
            f'bounds so corners will be interpolated and extrapolated from '
            f'cell centers instead.',
            stacklevel=2,
        )

    return interp_extrap_corners_2d(lat), interp_extrap_corners_2d(lon)


def _get_cf_bounds(ds, var_name, shape):
    """
    Get the CF ``bounds`` of the given variable as a numpy array with the
    expected shape, or ``None`` if they are not available
    """
    bounds_name = ds[var_name].attrs.get('bounds')
    if bounds_name is None:
        return None

    if bounds_name not in ds:
        warnings.warn(
            f'{var_name} has a CF bounds attribute "{bounds_name}" but no '
            f'such variable is present in the dataset.',
            stacklevel=3,
        )
        return None

    bounds = np.array(ds[bounds_name].values, float)
    if bounds.shape != shape:
        warnings.warn(
            f'The CF bounds variable {bounds_name} has shape '
            f'{bounds.shape}, not the expected {shape}.',
            stacklevel=3,
        )
        return None

    return bounds


def _bounds_tolerance(bounds):
    """A tolerance for comparing bounds, based on the size of the cells"""
    center = np.mean(bounds, axis=-1, keepdims=True)
    scale = np.max(np.abs(bounds - center))
    return 1e-6 * scale


def _corners_from_bounds_1d(bounds):
    """
    Convert CF bounds with shape ``(n, 2)`` to an array of ``n + 1`` corners,
    or return ``None`` if the bounds are not contiguous (in which case they
    cannot be described by a 1D array of corners)
    """
    tol = _bounds_tolerance(bounds)
    # bounds may be given in the direction of the coordinate or always from
    # the lower to the upper edge, so try both orders
    for flipped in [bounds, bounds[:, ::-1]]:
        if np.all(np.abs(flipped[:-1, 1] - flipped[1:, 0]) <= tol):
            return np.append(flipped[:, 0], flipped[-1, 1])

    return None


def _corners_from_bounds_2d(lat_bounds, lon_bounds):
    """
    Convert CF bounds with shape ``(ny, nx, 4)`` to arrays of corners with
    shape ``(ny + 1, nx + 1)``, or return ``None`` if neighboring cells do not
    share vertices (in which case they cannot be described by 2D arrays of
    corners)
    """
    # CF requires the vertices of each cell to be traversed in order around
    # the cell but does not say which vertex comes first or which direction
    # the traversal goes, so try each of the 8 possibilities.  Each candidate
    # gives the vertex indices of the lower-left, lower-right, upper-right and
    # upper-left corners of the cell in index space.  Requiring neighboring
    # cells to share vertices picks out the right candidate except on grids
    # too small to have neighbors in one or both directions, where the CF
    # recommendation (anticlockwise starting from the lower left) is assumed.
    candidates = []
    for base in ([0, 1, 2, 3], [0, 3, 2, 1]):
        for shift in range(4):
            candidates.append(base[shift:] + base[:shift])

    tol = max(_bounds_tolerance(lat_bounds), _bounds_tolerance(lon_bounds))

    for candidate in candidates:
        if not all(
            _vertices_are_shared(bounds, candidate, tol)
            for bounds in [lat_bounds, lon_bounds]
        ):
            continue

        lower_left, lower_right, upper_right, upper_left = candidate
        corners = []
        for bounds in [lat_bounds, lon_bounds]:
            ny, nx = bounds.shape[0], bounds.shape[1]
            corner = np.zeros((ny + 1, nx + 1))
            corner[:-1, :-1] = bounds[:, :, lower_left]
            corner[:-1, -1] = bounds[:, -1, lower_right]
            corner[-1, :-1] = bounds[-1, :, upper_left]
            corner[-1, -1] = bounds[-1, -1, upper_right]
            corners.append(corner)

        return corners[0], corners[1]

    return None


def _vertices_are_shared(bounds, candidate, tol):
    """
    Whether neighboring cells share vertices if the 4 vertices of each cell
    are in the order given by ``candidate`` (lower-left, lower-right,
    upper-right and upper-left in index space)
    """
    lower_left, lower_right, upper_right, upper_left = candidate
    # neighbors in x share their left and right vertices while neighbors in y
    # share their lower and upper vertices
    shared = [
        (bounds[:, :-1, lower_right], bounds[:, 1:, lower_left]),
        (bounds[:, :-1, upper_right], bounds[:, 1:, upper_left]),
        (bounds[:-1, :, upper_left], bounds[1:, :, lower_left]),
        (bounds[:-1, :, upper_right], bounds[1:, :, lower_right]),
    ]
    return all(
        np.all(np.abs(first - second) <= tol) for first, second in shared
    )


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

    attrs = ds.grid_corner_lat.attrs
    ds['grid_corner_lat'] = (('grid_size', 'grid_corners'), grid_corner_lat)
    ds.grid_corner_lat.attrs = attrs
    attrs = ds.grid_corner_lon.attrs
    ds['grid_corner_lon'] = (('grid_size', 'grid_corners'), grid_corner_lon)
    ds.grid_corner_lon.attrs = attrs


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
