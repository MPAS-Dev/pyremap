#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
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

import xarray

from pyremap import MpasCellMeshDescriptor, Remapper, get_polar_descriptor


# replace with the MPAS mesh name
inGridName = 'oQU240'

# replace with the path to the desired mesh or restart file
# As an example, use:
# https://web.lcrc.anl.gov/public/e3sm/inputdata/ocn/mpas-o/oQU240/ocean.QU.240km.151209.nc
inGridFileName = 'ocean.QU.240km.151209.nc'

inDescriptor = MpasCellMeshDescriptor(inGridFileName, inGridName)

# modify the size and resolution of the Antarctic grid as desired
outDescriptor = get_polar_descriptor(Lx=6000., Ly=5000., dx=10., dy=10.,
                                     projection='antarctic')
outGridName = outDescriptor.meshName

mappingFileName = f'map_{inGridName}_to_{outGridName}_bilinear.nc'

remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

# conservative remapping with 4 MPI tasks (using mpirun)
remapper.esmf_build_map(method='bilinear', mpi_tasks=4)

# select the SST at the initial time as an example data set
srcFileName = f'temp_{inGridName}.nc'
ds = xarray.open_dataset(inGridFileName)
dsOut = xarray.Dataset()
dsOut['temperature'] = ds['temperature'].isel(nVertLevels=0, Time=0)
dsOut.to_netcdf(srcFileName)

# do remapping with ncremap
outFileName = f'temp_{outGridName}_file.nc'
remapper.remap_file(srcFileName, outFileName)

# do remapping again, this time with python remapping
outFileName = f'temp_{outGridName}_array.nc'
dsOut = remapper.remap(dsOut)
dsOut.to_netcdf(outFileName)
