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
import netCDF4


def write_netcdf(ds, filename, fill_values=netCDF4.default_fillvals,
                 format='NETCDF4'):
    '''Write an xarray Dataset with NetCDF4 fill values where needed'''
    encoding_dict = {}
    varnames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for varname in varnames:
        is_numeric = numpy.issubdtype(ds[varname].dtype, numpy.number)
        if is_numeric and numpy.any(numpy.isnan(ds[varname])):
            dtype = ds[varname].dtype
            for fill_type in fill_values:
                if dtype == numpy.dtype(fill_type):
                    encoding_dict[varname] = \
                        {'_FillValue': fill_values[fill_type]}
                    break
        else:
            encoding_dict[varname] = {'_FillValue': None}

    ds.to_netcdf(filename, encoding=encoding_dict, format=format)
