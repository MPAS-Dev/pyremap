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
import pyproj
import xarray
from collections import OrderedDict

from pyremap.mesh_descriptor import MeshDescriptor


class ProjectionGridDescriptor(MeshDescriptor):
    """
    A class for describing a general logically rectangular grid that can be
    defined by a `pyproj` projection.

    Attributes
    ----------
    projection : pyproj.Proj
        The projection used to map from grid x-y space to latitude and
        longitude or ``None`` if the grid is a lat/lon grid

    lat_lon_projection : pyproj.Proj
        The reference projection used to get lat/lon from x/y
    """

    def __init__(self, projection=None, ds=None, filename=None, meshname=None,
                 x=None, y=None, xbnds=None, ybnds=None, xvarname='x',
                 yvarname='y', xdimname=None, ydimname=None, xunits=None,
                 yunits=None):
        """
        Constructor stores the projection

        Parameters
        ----------
        projection : pyproj.Proj, optional
            The projection used to map from grid x-y space to latitude and
            longitude.  If none is provided, the grid is a lat/lon grid

        ds : xarray.Dataset, optional
            A data set to read coordinates from.  One of ``ds`` or ``filename``
            must be provided

        filename : str, optional
            The path of the file containing the grid data.  One of ``ds`` or
            ``filename`` must be provided

        meshname : str, optional
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``).  If not
            provided, the data set in ``filename`` must have a global
            attribute ``meshname`` that will be used instead.

        meshname : str
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``)

        x, y : numpy array, optional
            One dimensional arrays defining the x and y coordinates of grid
            cell centers.

        xbnds, ybnds : numpy.arrays, optional
            2D arrays of bounds for each x and y cell, with size len(x) by 2
            and len(y) by 2, respectively.  These are interpolated and
            extrapolated from ``x`` and ``y`` if not provided.

        xvarname, yvarname : str, optional
            The name of the x and y variables in the grid file

        xdimname, ydimname : str, optional
            The name of the x and y dimensions, the same as xvarname and
            yvarname by default

        xunits, yunits : str, optional
            Will be taken from the ``units`` attribute of ``xvarname`` and
            ``yvarname`` if not explicitly provided, or 'm' if x and y
            and/or xbnds and ybnds are provided
        """

        super().__init__()

        self.projection = projection
        self.lat_lon_projection = pyproj.Proj(proj='latlong', datum='WGS84')

        coords_provided = ((x is not None and y is not None) or
                           (xbnds is not None and ybnds is not None))

        if ds is None and filename is None and not coords_provided:
            raise ValueError('Must provide ds, filename, x and y, or xbnds '
                             'and ybnds.')

        if coords_provided:
            if x is None:
                x = 0.5 * (xbnds[:, 0] + xbnds[:, 1])

            if y is None:
                y = 0.5 * (ybnds[:, 0] + ybnds[:, 1])

            if xbnds is None:
                xbnds = MeshDescriptor.interp_extrap_corner(x)

            if ybnds is None:
                ybnds = MeshDescriptor.interp_extrap_corner(y)

            if xdimname is None:
                xdimname = xvarname

            if ydimname is None:
                ydimname = yvarname

            if xunits is None:
                xunits = 'm'
            if yunits is None:
                yunits = 'm'

        else:
            if ds is None:
                ds = xarray.open_dataset(filename)

            if meshname is None and 'meshname' in ds.attrs:
                meshname = ds.attrs['meshname']

            # Get info from input file
            x = numpy.array(ds[xvarname].values, float)
            y = numpy.array(ds[yvarname].values, float)

            if 'bounds' in ds[xvarname].attrs:
                xbvarname = ds[xvarname].attrs['bounds']
                xbnds = numpy.array(ds[xbvarname].values, float)
            else:
                xbnds = MeshDescriptor.interp_extrap_corner(x)

            if 'bounds' in ds[yvarname].attrs:
                ybvarname = ds[yvarname].attrs['bounds']
                ybnds = numpy.array(ds[ybvarname].values, float)
            else:
                ybnds = MeshDescriptor.interp_extrap_corner(y)

            if xunits is None:
                if 'units' not in ds[xvarname].attrs:
                    raise ValueError('No units were provided for x.')
                xunits = ds[xvarname].attrs['units']
            if yunits is None:
                if 'units' not in ds[yvarname].attrs:
                    raise ValueError('No units were provided for y.')
                yunits = ds[yvarname].attrs['units']

            xdimname = ds[xvarname].dims[0]
            ydimname = ds[yvarname].dims[0]

        self._set_coords(x, y, xbnds, ybnds, xvarname, yvarname,
                         xdimname, ydimname, xunits, yunits, meshname)

        # Update history attribute of netCDF file
        if ds is not None and 'history' in ds.attrs:
            self.history = '\n'.join([ds.attrs['history'], self.history])

    def project_to_lon_lat(self, X, Y, xunits, yunits):
        """
        Given X and Y locations of points in a projection, returns the
        corresponding latitude and longitude of each point.

        Parameters
        ----------
        X, Y : numpy.array
            X and y arrays of points in in the projection

        xunits, yunits : str
            The units of X and Y if this is a lat/lon grid

        Returns
        -------
        Lat, Lon : numpy.array with same shape as X and Y
            the latitude and longitude in degrees of the points
        """

        if self.projection is None:
            if 'degree' in xunits:
                Lon = X
            else:
                Lon = numpy.rad2deg(X)
            if 'degree' in yunits:
                Lat = Y
            else:
                Lat = numpy.rad2deg(Y)
        else:
            Lon, Lat = pyproj.transform(self.projection,
                                        self.lat_lon_projection,
                                        X, Y)

        return Lon, Lat

    def _set_coords(self, x, y, xbnds, ybnds, xvarname, yvarname, xdimname,
                    ydimname, xunits, yunits, meshname):
        """
        Set up a coords dict with x, y, lat and lon
        """

        if meshname is None:
            dx = MeshDescriptor._round_res(abs(x[1] - x[0]))
            dy = MeshDescriptor._round_res(abs(y[1] - y[0]))
            if 'degree' in xunits and 'degree' in yunits:
                meshname = '{}x{}{}'.format(dx, dy, 'degrees')
            elif xunits == yunits:
                meshname = '{}x{}{}'.format(dx, dy, xunits)
            else:
                meshname = '{}{}x{}{}'.format(dx, xunits, dy, yunits)

        self.meshname = meshname

        if self.projection is None:
            # lat/lon grid
            lon_range = xbnds[-1, 1] - xbnds[0, 0]
            lat_range = ybnds[-1, 1] - ybnds[0, 0]
            if 'rad' in xunits:
                lon_range = numpy.rad2deg(lon_range)
            if 'rad' in yunits:
                lat_range = numpy.rad2deg(lat_range)

            # does the grid cover the full sphere?
            self.regional = (numpy.abs(lon_range - 360.) > 1e-10 or
                             numpy.abs(lat_range - 180.) > 1e-10)
        else:
            self.regional = True

        (X, Y) = numpy.meshgrid(x, y)
        (Lon, Lat) = self.project_to_lon_lat(X, Y, xunits, yunits)
        self.set_lon_lat_cell_centers(Lon, Lat, degrees=True)

        xbndsname = '{}_bnds'.format(xvarname)
        ybndsname = '{}_bnds'.format(yvarname)

        # make bounds into a 2D array as expected by set_lat_lon_vertices
        xb = numpy.zeros(len(x)+1)
        xb[0:-1] = xbnds[:, 0]
        xb[-1] = xbnds[-1, 1]

        yb = numpy.zeros(len(y)+1)
        yb[0:-1] = ybnds[:, 0]
        yb[-1] = ybnds[-1, 1]

        Xc, Yc = numpy.meshgrid(xb, yb)
        (Lonc, Latc) = self.project_to_lon_lat(Xc, Yc, xunits, yunits)

        self.set_lon_lat_vertices(Lonc, Latc, degrees=True)

        if self.projection is None:
            xlong = 'longitude'
            xstd = 'longitude'
            ylong = 'latitude'
            ystd = 'latitude'
        else:
            xlong = 'x coordinate of projection'
            xstd = 'projection_x_coordinate'
            ylong = 'y coordinate of projection'
            ystd = 'projection_y_coordinate'

        self.coords = OrderedDict(
            [(xvarname, {'dims': (xdimname,),
                         'data': x,
                         'attrs': {'units': xunits,
                                   'long_name': xlong,
                                   'standard_name': xstd,
                                   'bounds': xbndsname,
                                   'axis': 'X'}}),
             (xbndsname, {'dims': (xdimname, 'nbnds'),
                          'data': xbnds}),
             (yvarname, {'dims': (ydimname,),
                         'data': y,
                         'attrs': {'units': yunits,
                                   'long_name': ylong,
                                   'standard_name': ystd,
                                   'bounds': ybndsname,
                                   'axis': 'Y'}}),
             (ybndsname, {'dims': (ydimname, 'nbnds'),
                          'data': ybnds})])
        if self.projection is None:
            self.lon_lat_coords = [yvarname, xvarname]
        else:
            self.lon_lat_coords = ['lat', 'lon']
            self.coords['lat'] = {'dims': (ydimname, xdimname),
                                  'data': Lat,
                                  'attrs': {'units': 'degrees_north',
                                            'long_name': 'latitude',
                                            'standard_name': 'latitude',
                                            'bounds': 'lat_vertices'}}
            self.coords['lon'] = {'dims': (ydimname, xdimname),
                                  'data': Lon,
                                  'attrs': {'units': 'degrees_east',
                                            'long_name': 'longitude',
                                            'standard_name': 'longitude',
                                            'bounds': 'lon_vertices'}}
            self.coords['lat_vertices'] = \
                {'dims': (ydimname, xdimname, 'nv'),
                 'data': MeshDescriptor._unwrap_corners(Latc)}
            self.coords['lon_vertices'] = \
                {'dims': (ydimname, xdimname, 'nv'),
                 'data': MeshDescriptor._unwrap_corners(Lonc)}

        # y first, then x
        self.sizes = OrderedDict([(ydimname, len(y)), (xdimname, len(x))])
