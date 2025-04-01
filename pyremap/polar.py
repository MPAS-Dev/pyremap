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
import numpy
import pyproj
import xarray

from pyremap.descriptor import ProjectionGridDescriptor


def get_arctic_stereographic_projection():
    """
    Get a projection for an Arctic stereographic comparison grid

    Returns
    -------
    projection : ``pyproj.Proj`` object
        The projection
    """
    # Authors
    # -------
    # Milena Veneziani

    projection = pyproj.Proj(
        '+proj=stere +lat_ts=75.0 +lat_0=90 +lon_0=0.0 '
        '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'
    )

    return projection


def get_antarctic_stereographic_projection():
    """
    Get a projection for an Antarctic steregraphic grid
    """

    projection = pyproj.Proj(
        '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
        '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'
    )

    return projection


def get_polar_descriptor_from_file(filename, projection='antarctic'):
    """
    Get a descriptor of a polar stereographic grid used for remapping

    Parameters
    ----------
    filename :  str
        A file containing x and y coordinates for the grid

    projection : {'arctic', 'antarctic', pyproj.Proj}
        The projection to use.  'arctic' and 'antarctic' are polar
        stereographic projections with reference latitude at +/- 71 degrees

    Returns
    -------
    descriptor : ``ProjectionGridDescriptor`` object
        A descriptor of the polar grid
    """
    ds_in = xarray.open_dataset(filename)
    x = ds_in.x.values
    y = ds_in.y.values
    dx = int((x[1] - x[0]) / 1000.0)
    lx = int((x[-1] - x[0]) / 1000.0)
    ly = int((y[-1] - y[0]) / 1000.0)

    mesh_name = f'{lx}x{ly}km_{dx}km_antarctic_stereo'

    projection = _get_projection(projection)

    descriptor = ProjectionGridDescriptor.create(projection, x, y, mesh_name)

    return descriptor


def get_polar_descriptor(lx, ly, dx, dy, projection='antarctic'):
    """
    Get a descriptor of a polar stereographic grid used for remapping

    Parameters
    ----------
    lx, ly :  double
        Size of the domain in x and y in km

    dx, dy : double
        Resolution of the grid in km

    projection : {'arctic', 'antarctic', pyproj.Proj}
        The projection to use.  'arctic' and 'antarctic' are polar
        stereographic projections with reference latitude at +/- 71 degrees

    Returns
    -------
    descriptor : ``ProjectionGridDescriptor`` object
        A descriptor of the Antarctic grid
    """

    upper_proj = projection[0].upper() + projection[1:]

    mesh_name = f'{lx}x{ly}km_{dx}km_{upper_proj}_stereo'

    x_max = 0.5 * lx * 1e3
    nx = int(lx / dx) + 1
    x = numpy.linspace(-x_max, x_max, nx)

    y_max = 0.5 * ly * 1e3
    ny = int(ly / dy) + 1
    y = numpy.linspace(-y_max, y_max, ny)

    projection = _get_projection(projection)

    descriptor = ProjectionGridDescriptor.create(projection, x, y, mesh_name)

    return descriptor


def to_polar(points):
    projection = get_antarctic_stereographic_projection()
    lat_lon_projection = pyproj.Proj(proj='latlong', datum='WGS84')

    transformer = pyproj.Transformer.from_proj(lat_lon_projection, projection)
    x, y = transformer.transform(points[:, 0], points[:, 1], radians=False)
    points[:, 0] = x
    points[:, 1] = y
    return points


def from_polar(points):
    projection = get_antarctic_stereographic_projection()
    lat_lon_projection = pyproj.Proj(proj='latlong', datum='WGS84')

    transformer = pyproj.Transformer.from_proj(projection, lat_lon_projection)
    lon, lat = transformer.transform(points[:, 0], points[:, 1], radians=False)
    points[:, 0] = lon
    points[:, 1] = lat
    return points


def _get_projection(projection):
    if isinstance(projection, str):
        if projection == 'arctic':
            projection = get_arctic_stereographic_projection()
        elif projection == 'antarctic':
            projection = get_antarctic_stereographic_projection()
        else:
            raise ValueError(f'Bad projection name {projection}')
    return projection
