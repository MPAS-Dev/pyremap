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

from pyremap.projection_grid_descriptor import ProjectionGridDescriptor


class LonLatGridDescriptor(ProjectionGridDescriptor):
    """
    A class for describing a latitude/longitude tensor grid
    """

    def __init__(self, ds=None, filename=None, meshname=None, lon=None,
                 lat=None, lonbnds=None, latbnds=None, lonvarname='lon',
                 latvarname='lat', londimname=None, latdimname=None,
                 degrees=None):
        """
        Constructor stores the projection

        Parameters
        ----------
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

        lon, lat : numpy array, optional
            One dimensional arrays defining the lon and lat coordinates of grid
            cell centers.

        lonbnds, latbnds : numpy.arrays, optional
            2D arrays of bounds for each lon and lat cell, with size len(lon)
            by 2 and len(lat) by 2, respectively.  These are interpolated and
            extrapolated from ``lon`` and ``lat`` if not provided.

        lonvarname, latvarname : str, optional
            The name of the lon and lat variables in the grid file

        londimname, latdimname : str, optional
            The name of the lon and lat dimensions, the same as lonvarname and
            latvarname by default

        degrees : bool, optional
            If the data is in degrees (vs. radians).  If not provided, it will
            be inferred from the ``units`` attribute on ``lonvarname``.
        """

        if degrees is None:
            xunits = None
            yunits = None
        elif degrees:
            xunits = 'degrees_east'
            yunits = 'degrees_north'
        else:
            xunits = 'radians_east'
            yunits = 'radians_north'

        super().__init__(projection=None, ds=ds, filename=filename,
                         meshname=meshname, x=lon, y=lat, xbnds=lonbnds,
                         ybnds=latbnds, xvarname=lonvarname,
                         yvarname=latvarname, xdimname=londimname,
                         ydimname=latdimname,  xunits=xunits, yunits=yunits)

    @classmethod
    def create_constant_spacing(cls, dlon, dlat, lon_min=-180., lon_max=180.,
                                lat_min=-90., lat_max=90.):
        """
        Create a descriptor of a lat/lon grid with constant grid spacing

        Parameters
        ----------
        dlon, dlat : float
            Resolution of the lat/lon grid in degrees

        lon_min, lon_max, lat_min, lat_max : float, optional
            Bounds of the lat/lon grid in degrees

        Returns
        -------
        descriptor : ``LatLonGridDescriptor`` object
            A descriptor of the lat/lon grid
        """
        nlat = int((lat_max - lat_min) / dlat)
        nlon = int((lon_max - lon_min) / dlon)
        latbnds = numpy.zeros((nlat, 2))
        latbnds[:, 0] = numpy.linspace(lat_min, lat_max-dlat, nlat)
        latbnds[:, 1] = numpy.linspace(lat_min+dlat, lat_max, nlat)
        lonbnds = numpy.zeros((nlon, 2))
        lonbnds[:, 0] = numpy.linspace(lon_min, lon_max-dlon, nlon)
        lonbnds[:, 1] = numpy.linspace(lon_min+dlon, lon_max, nlon)
        lat = 0.5*(latbnds[:, 0] + latbnds[:, 1])
        lon = 0.5*(lonbnds[:, 0] + lonbnds[:, 1])

        descriptor = LonLatGridDescriptor(
            lon=lon, lat=lat, lonbnds=lonbnds, latbnds=latbnds, degrees=True)

        return descriptor
