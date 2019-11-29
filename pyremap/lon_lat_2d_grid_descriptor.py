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

from pyremap.mesh_descriptor import MeshDescriptor


class LonLat2DGridDescriptor(MeshDescriptor):
    """
    A class for describing a lat-lon grid that may not be a tensor grid
    (lat/lon are 2D arrays).  The grid is assumed to be regional, since this
    is difficult to determine just from the lat/lon values.  The calling code
    should set ``regional = False`` for global grids with 2D lat/lon
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, ds=None, filename=None, meshname=None, lonvarname='lon',
                 latvarname='lat', lonvertname=None, latvertname=None,
                 units=None):
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

        lonvarname, latvarname : str, optional
            The name of the lon and lat variables for cell centers from the grid
             file

        lonvertname, latvertname : str, optional
            The name of the lon and lat variables for cell vertices in the grid
            file.  These must be one element larger in each dimension than
            the variables ``lonvarname`` and ``latvarname`` if provided.  If
            not provided, they are interpolated and extrapolated from
            ``lonvarname`` and ``latvarname``.

        units : str, optional
            Will be taken from the ``units`` attribute of ``lonvarname`` if not
            explicitly provided
        """

        super().__init__()

        if ds is None:
            ds = xarray.open_dataset(filename)

        if meshname is None and 'meshname' in ds.attrs:
            meshname = ds.attrs['meshname']

        # Get info from input file
        lon = numpy.array(ds[lonvarname].values, float)
        lat = numpy.array(ds[latvarname].values, float)
        if units is None:
            if 'degree' in ds[latvarname].units:
                units = 'degrees'
            elif 'rad' in ds[latvarname].units:
                units = 'radians'
            else:
                raise ValueError("Couldn't figure out units {}".format(
                    ds[latvarname].units))

        if lonvertname is None:
            lonvert = MeshDescriptor.interp_extrap_corners_2d(lon)
        else:
            lonvert = numpy.array(ds[lonvertname].values, float)

        if latvertname is None:
            latvert = MeshDescriptor.interp_extrap_corners_2d(lat)
        else:
            latvert = numpy.array(ds[latvertname].values, float)

        self._set_coords(lon, lat, lonvert, latvert, latvarname, lonvarname,
                         ds[latvarname].dims[0], ds[latvarname].dims[1],
                         units, meshname)

        # Update history attribute of netCDF file
        if 'history' in ds.attrs:
            self.history = '\n'.join([ds.attrs['history'], self.history])

    def _set_coords(self, lon, lat, lonvert, latvert, lonvarname,
                    latvarname, xdimname,  ydimname, units, meshname):
        """
        Set up a coords dict with x, y, lat and lon
        """
        dlon = MeshDescriptor._round_res(abs(lon[0, 1] - lon[0, 0]))
        dlat = MeshDescriptor._round_res(abs(lat[1, 0] - lat[0, 0]))

        if meshname is None:
            meshname = '{}x{}{}'.format(dlon, dlat, units)

        self.meshname = meshname

        self.regional = True

        self.set_lon_lat_cell_centers(lon, lat, degrees=True)

        self.set_lon_lat_vertices(lonvert, latvert)

        self.coords = {latvarname: {'dims': (ydimname, xdimname),
                                    'data': lat,
                                    'attrs': {'units': units}},
                       lonvarname: {'dims': (ydimname, xdimname),
                                    'data': lon,
                                    'attrs': {'units': units}}}

        self.sizes = {xdimname: lat.shape[1], ydimname: lat.shape[0]}
