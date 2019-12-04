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
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

"""
Creates a mapping file that can be used with ncremap (NCO) to remap MPAS files
to a latitude/longitude grid.

Usage: Copy this script into the main MPAS-Analysis directory (up one level).
Modify the grid name, the path to the MPAS grid file and the output grid
resolution.
"""

import xarray

from pyremap import MpasMeshDescriptor, Remapper, LonLatGridDescriptor, \
    write_netcdf

# replace with the MPAS mesh name
inGridName = 'oQU240'

# replace with the path to the desired mesh or restart file
# As an example, use:
# https://web.lcrc.anl.gov/public/e3sm/inputdata/ocn/mpas-o/oQU240/ocean.QU.240km.151209.nc
inGridFileName = 'ocean.QU.240km.151209.nc'

inDescriptor = MpasMeshDescriptor(inGridFileName, inGridName)

# modify the resolution of the global lat-lon grid as desired
outDescriptor = LonLatGridDescriptor.create_constant_spacing(dlon=0.1, dlat=0.1)
outGridName = outDescriptor.meshname

mappingFileName = 'map_{}_to_{}_bilinear.nc'.format(inGridName, outGridName)


remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

remapper.build_mapping_file(method='bilinear', mpitasks=8)

outFileName = 'sst_{}.nc'.format(outGridName)
ds = xarray.open_dataset(inGridFileName)
dsOut = xarray.Dataset()
mask = ds.maxLevelCell > 0
dsOut['sst'] = ds['temperature'].isel(nVertLevels=0).where(mask)
dsOut = remapper.remap(dsOut, renorm_thresh=0.01)
write_netcdf(dsOut, outFileName, always_fill=False)
