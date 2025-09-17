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
Creates a mapping file that can be used with ncremap (NCO) to remap MPAS files
to a latitude/longitude grid.

Usage: Copy this script into the main MPAS-Analysis directory (up one level).
Modify the grid name, the path to the MPAS grid file and the output grid
resolution.
"""

import xarray as xr

from pyremap import Remapper, get_polar_descriptor

# replace with the MPAS mesh name
src_mesh_name = 'oQU240'

# replace with the path to the desired mesh or restart file
# As an example, use:
# https://web.lcrc.anl.gov/public/e3sm/inputdata/ocn/mpas-o/oQU240/ocean.QU.240km.151209.nc
src_mesh_filename = 'ocean.QU.240km.151209.nc'

# Configure the remapper to build the mapping file using 4 MPI tasks and
# bilinear remapping
remapper = Remapper(ntasks=4, method='bilinear')

remapper.src_from_mpas(filename=src_mesh_filename, mesh_name=src_mesh_name)

# modify the size and resolution of the Antarctic grid as desired
remapper.dst_descriptor = get_polar_descriptor(
    lx=6000.0, ly=5000.0, dx=10.0, dy=10.0, projection='antarctic'
)

dst_grid_name = remapper.dst_descriptor.mesh_name

# build the mapping file
remapper.build_map()

# select the SST at the initial time as an example data set
src_filename = f'temp_{src_mesh_name}.nc'
ds = xr.open_dataset(src_mesh_filename)
ds_out = xr.Dataset()
ds_out['temperature'] = ds['temperature'].isel(nVertLevels=0, Time=0)
ds_out.to_netcdf(src_filename)

# do remapping with ncremap
out_filename = f'temp_{dst_grid_name}_file.nc'
remapper.ncremap(src_filename, out_filename)

# do remapping again, this time with python remapping
out_filename = f'temp_{dst_grid_name}_array.nc'
ds_out = remapper.remap_numpy(ds_out)
ds_out.to_netcdf(out_filename)
