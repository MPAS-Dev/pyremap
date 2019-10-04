import pyproj
import xarray
import numpy

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

    projection = pyproj.Proj('+proj=stere +lat_ts=75.0 +lat_0=90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

    return projection


def get_antarctic_stereographic_projection():
    """
    Get a projection for an Antarctic steregraphic grid
    """

    projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

    return projection


def get_polar_descriptor_from_file(fileName, projection='antarctic'):
    """
    Get a descriptor of a polar stereographic grid used for remapping

    Parameters
    ----------
    fileName :  str
        A file containing x and y coordinates for the grid

    projection : {'arctic', 'antarctic', pyproj.Proj}
        The projection to use.  'arctic' and 'antarctic' are polar
        stereographic projections with reference latitude at +/- 71 degrees

    Returns
    -------
    descriptor : ``ProjectionGridDescriptor`` object
        A descriptor of the polar grid
    """
    dsIn = xarray.open_dataset(fileName)
    x = dsIn.x.values
    y = dsIn.y.values
    dx = int((x[1]-x[0])/1000.)
    Lx = int((x[-1] - x[0])/1000.)
    Ly = int((y[-1] - y[0])/1000.)

    meshName = '{}x{}km_{}km_Antarctic_stereo'.format(Lx, Ly, dx)

    projection = _get_projection(projection)

    descriptor = ProjectionGridDescriptor.create(projection, x, y, meshName)

    return descriptor


def get_polar_descriptor(Lx, Ly, dx, dy, projection='antarctic'):
    """
    Get a descriptor of a polar stereographic grid used for remapping

    Parameters
    ----------
    Lx, Ly :  double
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

    meshName = '{}x{}km_{}km_Antarctic_stereo'.format(Lx, Ly, dx)

    xMax = 0.5 * Lx * 1e3
    nx = int(Lx / dx) + 1
    x = numpy.linspace(-xMax, xMax, nx)

    yMax = 0.5 * Ly * 1e3
    ny = int(Ly / dy) + 1
    y = numpy.linspace(-yMax, yMax, ny)

    projection = _get_projection(projection)

    descriptor = ProjectionGridDescriptor.create(projection, x, y, meshName)

    return descriptor


def to_polar(points):

    projection = get_antarctic_stereographic_projection()
    latLonProjection = pyproj.Proj(proj='latlong', datum='WGS84')

    x, y = pyproj.transform(latLonProjection, projection, points[:, 0],
                            points[:, 1], radians=False)
    points[:, 0] = x
    points[:, 1] = y
    return points


def from_polar(points):

    projection = get_antarctic_stereographic_projection()
    latLonProjection = pyproj.Proj(proj='latlong', datum='WGS84')

    lon, lat = pyproj.transform(projection, latLonProjection, points[:, 0],
                                points[:, 1], radians=False)
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
            raise ValueError('Bad projection name {}'.format(projection))
    return projection
