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
"""
Unit tests for using CF ``bounds`` to locate grid-cell corners.
"""

import numpy as np
import pyproj
import pytest
import xarray as xr

from pyremap import (
    LatLon2DGridDescriptor,
    LatLonGridDescriptor,
    ProjectionGridDescriptor,
)
from pyremap.descriptor import (
    get_corners_1d,
    get_corners_2d,
    interp_extrap_corner,
    interp_extrap_corners_2d,
)

# corners of a grid with cells that vary in size, so that the corners are not
# the same as those interpolated and extrapolated from the cell centers
LAT_CORNER = np.array([-90.0, -60.0, -10.0, 20.0, 30.0, 90.0])
LON_CORNER = np.array([-180.0, -100.0, -30.0, 0.0, 45.0, 90.0, 180.0])


def _centers(corner):
    return 0.5 * (corner[:-1] + corner[1:])


def _bounds_1d(corner):
    """CF bounds of shape (n, 2) in the direction of the coordinate"""
    return np.stack((corner[:-1], corner[1:]), axis=-1)


def _lat_lon_dataset(lat_bounds=None, lon_bounds=None):
    """A 1D lat-lon dataset, optionally with CF bounds"""
    lat = _centers(LAT_CORNER)
    lon = _centers(LON_CORNER)
    ds = xr.Dataset(
        coords={
            'lat': ('lat', lat, {'units': 'degrees_north'}),
            'lon': ('lon', lon, {'units': 'degrees_east'}),
        }
    )
    if lat_bounds is not None:
        ds['lat_bnds'] = (('lat', 'nbnd'), lat_bounds)
        ds.lat.attrs['bounds'] = 'lat_bnds'
    if lon_bounds is not None:
        ds['lon_bnds'] = (('lon', 'nbnd'), lon_bounds)
        ds.lon.attrs['bounds'] = 'lon_bnds'
    return ds


def _lat_lon_2d_dataset(order=(0, 1, 2, 3), lat_corner=None, lon_corner=None):
    """
    A 2D lat-lon dataset with CF bounds, with the 4 vertices of each cell in
    the order given by ``order`` (lower-left, lower-right, upper-right and
    upper-left in index space)
    """
    if lat_corner is None:
        lat_corner = LAT_CORNER
    if lon_corner is None:
        lon_corner = LON_CORNER
    lon_corner_2d, lat_corner_2d = np.meshgrid(lon_corner, lat_corner)
    lat = 0.25 * (
        lat_corner_2d[:-1, :-1]
        + lat_corner_2d[:-1, 1:]
        + lat_corner_2d[1:, 1:]
        + lat_corner_2d[1:, :-1]
    )
    lon = 0.25 * (
        lon_corner_2d[:-1, :-1]
        + lon_corner_2d[:-1, 1:]
        + lon_corner_2d[1:, 1:]
        + lon_corner_2d[1:, :-1]
    )

    ds = xr.Dataset()
    for var_name, corner_2d, center in [
        ('lat2d', lat_corner_2d, lat),
        ('lon2d', lon_corner_2d, lon),
    ]:
        # lower-left, lower-right, upper-right, upper-left
        vertices = [
            corner_2d[:-1, :-1],
            corner_2d[:-1, 1:],
            corner_2d[1:, 1:],
            corner_2d[1:, :-1],
        ]
        bounds = np.zeros((center.shape[0], center.shape[1], 4))
        for vertex_index, corner_index in enumerate(order):
            bounds[:, :, vertex_index] = vertices[corner_index]

        units = 'degrees_north' if var_name == 'lat2d' else 'degrees_east'
        ds[var_name] = (('y', 'x'), center, {'units': units})
        ds[f'{var_name}_bnds'] = (('y', 'x', 'nv'), bounds)
        ds[var_name].attrs['bounds'] = f'{var_name}_bnds'

    return ds, lat_corner_2d, lon_corner_2d


def test_corners_1d_from_bounds():
    """CF bounds are used in place of interpolation and extrapolation"""
    ds = _lat_lon_dataset(
        lat_bounds=_bounds_1d(LAT_CORNER), lon_bounds=_bounds_1d(LON_CORNER)
    )

    np.testing.assert_allclose(get_corners_1d(ds, 'lat'), LAT_CORNER)
    np.testing.assert_allclose(get_corners_1d(ds, 'lon'), LON_CORNER)

    # the point of the test data is that these are not the same as the
    # interpolated and extrapolated corners
    assert not np.allclose(
        interp_extrap_corner(ds.lat.values), LAT_CORNER, atol=1e-10
    )


def test_corners_1d_no_bounds():
    """Without CF bounds, corners are interpolated and extrapolated"""
    ds = _lat_lon_dataset()

    np.testing.assert_allclose(
        get_corners_1d(ds, 'lat'), interp_extrap_corner(ds.lat.values)
    )


def test_corners_1d_descending():
    """CF bounds are used for coordinates that decrease with index"""
    lat_corner = LAT_CORNER[::-1]
    lat = _centers(lat_corner)
    ds = xr.Dataset(coords={'lat': ('lat', lat, {'units': 'degrees_north'})})
    ds['lat_bnds'] = (('lat', 'nbnd'), _bounds_1d(lat_corner))
    ds.lat.attrs['bounds'] = 'lat_bnds'

    np.testing.assert_allclose(get_corners_1d(ds, 'lat'), lat_corner)


def test_corners_1d_descending_min_max_bounds():
    """
    CF bounds are used for a descending coordinate whose bounds are always
    given from the lower to the upper edge
    """
    lat_corner = LAT_CORNER[::-1]
    lat = _centers(lat_corner)
    # each pair is [min, max] rather than in the direction of the coordinate
    bounds = _bounds_1d(lat_corner)[:, ::-1]
    ds = xr.Dataset(coords={'lat': ('lat', lat, {'units': 'degrees_north'})})
    ds['lat_bnds'] = (('lat', 'nbnd'), bounds)
    ds.lat.attrs['bounds'] = 'lat_bnds'

    np.testing.assert_allclose(get_corners_1d(ds, 'lat'), lat_corner)


def test_corners_1d_noncontiguous_bounds():
    """Bounds with gaps between cells can't be used, so we fall back"""
    bounds = _bounds_1d(LAT_CORNER)
    # shrink each cell so neighbors no longer share an edge
    center = np.mean(bounds, axis=-1, keepdims=True)
    bounds = center + 0.9 * (bounds - center)
    ds = _lat_lon_dataset(lat_bounds=bounds)

    with pytest.warns(UserWarning, match='not contiguous'):
        corner = get_corners_1d(ds, 'lat')

    np.testing.assert_allclose(corner, interp_extrap_corner(ds.lat.values))


def test_corners_1d_missing_bounds_variable():
    """A ``bounds`` attribute pointing at nothing is a warning, not an error"""
    ds = _lat_lon_dataset()
    ds.lat.attrs['bounds'] = 'lat_bnds'

    with pytest.warns(UserWarning, match='no such variable'):
        corner = get_corners_1d(ds, 'lat')

    np.testing.assert_allclose(corner, interp_extrap_corner(ds.lat.values))


def test_corners_1d_wrong_bounds_shape():
    """Bounds that aren't (n, 2) are a warning, not an error"""
    ds = _lat_lon_dataset()
    ds['lat_bnds'] = (('lat',), LAT_CORNER[:-1])
    ds.lat.attrs['bounds'] = 'lat_bnds'

    with pytest.warns(UserWarning, match='shape'):
        corner = get_corners_1d(ds, 'lat')

    np.testing.assert_allclose(corner, interp_extrap_corner(ds.lat.values))


@pytest.mark.parametrize(
    'order',
    [
        (0, 1, 2, 3),  # counterclockwise from the lower left
        (1, 2, 3, 0),  # counterclockwise from the lower right
        (0, 3, 2, 1),  # clockwise from the lower left
        (2, 1, 0, 3),  # clockwise from the upper right
    ],
)
def test_corners_2d_from_bounds(order):
    """
    2D CF bounds are used whatever order the vertices of each cell are
    traversed in
    """
    ds, lat_corner_2d, lon_corner_2d = _lat_lon_2d_dataset(order=order)

    lat_corner, lon_corner = get_corners_2d(ds, 'lat2d', 'lon2d')

    np.testing.assert_allclose(lat_corner, lat_corner_2d)
    np.testing.assert_allclose(lon_corner, lon_corner_2d)

    # again, these should not match what interpolation and extrapolation gives
    assert not np.allclose(
        interp_extrap_corners_2d(ds.lat2d.values), lat_corner_2d, atol=1e-10
    )


def test_corners_2d_no_bounds():
    """Without CF bounds, 2D corners are interpolated and extrapolated"""
    ds, _, _ = _lat_lon_2d_dataset()
    ds = ds.drop_vars(['lat2d_bnds', 'lon2d_bnds'])
    del ds.lat2d.attrs['bounds']
    del ds.lon2d.attrs['bounds']

    lat_corner, lon_corner = get_corners_2d(ds, 'lat2d', 'lon2d')

    np.testing.assert_allclose(
        lat_corner, interp_extrap_corners_2d(ds.lat2d.values)
    )
    np.testing.assert_allclose(
        lon_corner, interp_extrap_corners_2d(ds.lon2d.values)
    )


def test_corners_2d_unshared_vertices():
    """
    2D bounds that neighboring cells don't share can't be described by corner
    arrays, so we warn and fall back
    """
    ds, _, _ = _lat_lon_2d_dataset()
    bounds = ds.lat2d_bnds.values
    center = np.mean(bounds, axis=-1, keepdims=True)
    ds['lat2d_bnds'] = (ds.lat2d_bnds.dims, center + 0.9 * (bounds - center))

    with pytest.warns(UserWarning, match='do not share vertices'):
        lat_corner, lon_corner = get_corners_2d(ds, 'lat2d', 'lon2d')

    np.testing.assert_allclose(
        lat_corner, interp_extrap_corners_2d(ds.lat2d.values)
    )


def test_corners_2d_bounds_on_one_coord_only():
    """Bounds on only one of the two coordinates are not enough"""
    ds, _, _ = _lat_lon_2d_dataset()
    ds = ds.drop_vars(['lon2d_bnds'])
    del ds.lon2d.attrs['bounds']

    with pytest.warns(UserWarning, match='Only one of'):
        lat_corner, lon_corner = get_corners_2d(ds, 'lat2d', 'lon2d')

    np.testing.assert_allclose(
        lat_corner, interp_extrap_corners_2d(ds.lat2d.values)
    )


def test_lat_lon_descriptor_honors_bounds():
    """``LatLonGridDescriptor.read()`` uses CF bounds"""
    ds = _lat_lon_dataset(
        lat_bounds=_bounds_1d(LAT_CORNER), lon_bounds=_bounds_1d(LON_CORNER)
    )

    descriptor = LatLonGridDescriptor.read(ds=ds, mesh_name='test')

    np.testing.assert_allclose(descriptor.lat_corner, LAT_CORNER)
    np.testing.assert_allclose(descriptor.lon_corner, LON_CORNER)


def test_lat_lon_descriptor_scrip_from_bounds(tmp_path):
    """The corners in the SCRIP file come from the CF bounds"""
    ds = _lat_lon_dataset(
        lat_bounds=_bounds_1d(LAT_CORNER), lon_bounds=_bounds_1d(LON_CORNER)
    )

    descriptor = LatLonGridDescriptor.read(ds=ds, mesh_name='test')
    scrip_filename = str(tmp_path / 'scrip.nc')
    descriptor.to_scrip(scrip_filename)

    with xr.open_dataset(scrip_filename) as ds_scrip:
        corner_lat = ds_scrip.grid_corner_lat.values
        corner_lon = ds_scrip.grid_corner_lon.values

    nlat = len(LAT_CORNER) - 1
    nlon = len(LON_CORNER) - 1
    assert corner_lat.shape == (nlat * nlon, 4)
    # the first cell is bounded by the first two corners in each direction
    np.testing.assert_allclose(
        np.sort(np.unique(corner_lat[0, :])), LAT_CORNER[0:2]
    )
    np.testing.assert_allclose(
        np.sort(np.unique(corner_lon[0, :])), LON_CORNER[0:2]
    )


def test_lat_lon_2d_descriptor_honors_bounds():
    """``LatLon2DGridDescriptor.read()`` uses CF bounds"""
    ds, lat_corner_2d, lon_corner_2d = _lat_lon_2d_dataset()

    descriptor = LatLon2DGridDescriptor.read(
        ds=ds, lat_var_name='lat2d', lon_var_name='lon2d', mesh_name='test'
    )

    np.testing.assert_allclose(descriptor.lat_corner, lat_corner_2d)
    np.testing.assert_allclose(descriptor.lon_corner, lon_corner_2d)


def test_projection_descriptor_honors_bounds(tmp_path):
    """``ProjectionGridDescriptor.read()`` uses CF bounds on x and y"""
    x_corner = 1e3 * np.array([-300.0, -100.0, 0.0, 50.0, 300.0])
    y_corner = 1e3 * np.array([-200.0, -50.0, 100.0, 400.0])

    ds = xr.Dataset(
        coords={
            'x': ('x', _centers(x_corner), {'units': 'meters'}),
            'y': ('y', _centers(y_corner), {'units': 'meters'}),
        }
    )
    ds['x_bnds'] = (('x', 'nbnd'), _bounds_1d(x_corner))
    ds['y_bnds'] = (('y', 'nbnd'), _bounds_1d(y_corner))
    ds.x.attrs['bounds'] = 'x_bnds'
    ds.y.attrs['bounds'] = 'y_bnds'
    filename = str(tmp_path / 'proj_grid.nc')
    ds.to_netcdf(filename)

    projection = pyproj.Proj(
        '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 '
        '+x_0=0.0 +y_0=0.0 +ellps=WGS84'
    )
    descriptor = ProjectionGridDescriptor.read(
        projection=projection, filename=filename, mesh_name='test'
    )

    np.testing.assert_allclose(descriptor.x_corner, x_corner)
    np.testing.assert_allclose(descriptor.y_corner, y_corner)
