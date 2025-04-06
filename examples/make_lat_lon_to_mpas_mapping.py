#!/usr/bin/env python
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
Creates a mapping file that can be used with ncremap (NCO) to remap data from
a latitude/longitude grid to an MPAS-Ocean grid.

Usage: Modify the grid name, the path to the MPAS grid file, and the input grid
resolution.
"""

import numpy as np
import xarray as xr

from pyremap import Remapper, get_lat_lon_descriptor

# replace with the MPAS mesh name
dst_mesh_name = 'oQU240'

# replace with the path to the desired mesh or restart file
dst_mesh_filename = 'ocean.QU.240km.151209.nc'

# Configure the remapper to build the mapping file using 1 MPI task and
# bilinear remapping
remapper = Remapper(ntasks=1, method='bilinear')

# modify the resolution of the global lat-lon grid as desired
src_descriptor = get_lat_lon_descriptor(dlon=0.5, dlat=0.5)
remapper.src_descriptor = src_descriptor

remapper.dst_from_mpas(filename=dst_mesh_filename, mesh_name=dst_mesh_name)

src_grid_name = remapper.src_descriptor.mesh_name

# build the mapping file
remapper.build_map()

# create some made-up data on the lat-lon grid
lon2d, lat2d = np.meshgrid(src_descriptor.lon, src_descriptor.lat)
data = np.sin(np.radians(lat2d)) * np.cos(np.radians(lon2d))

ds = xr.Dataset()
ds['temperature'] = (('lat', 'lon'), data)
ds['lat'] = ('lat', src_descriptor.lat)
ds['lon'] = ('lon', src_descriptor.lon)

src_filename = f'temp_{src_grid_name}.nc'
ds.to_netcdf(src_filename)

# do remapping with ncremap
out_filename = f'temp_{dst_mesh_name}_file.nc'
remapper.ncremap(src_filename, out_filename)

# do remapping again, this time with python remapping
out_filename = f'temp_{dst_mesh_name}_array.nc'
ds_out = remapper.remap_numpy(ds)
ds_out.to_netcdf(out_filename)
