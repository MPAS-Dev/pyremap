#!/usr/bin/env python

'''
Creates a mapping file and remaps the variables from an input file on an
Antarctic grid to another grid with the same extent but a different resolution.
The mapping file can be used with ncremap (NCO) to remap MPAS files between
these same gridss.

Usage: Copy this script into the main MPAS-Analysis directory (up one level).
'''

import numpy
import xarray
import argparse

from pyremap import Remapper, ProjectionGridDescriptor
from pyremap.polar import get_antarctic_stereographic_projection


parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', dest='inFileName', required=True, type=str,
                    help="Input file name")
parser.add_argument('-o', dest='outFileName', required=True, type=str,
                    help="Output file name")
parser.add_argument('-r', dest='resolution', required=True, type=float,
                    help="Output resolution")
args = parser.parse_args()


dsIn = xarray.open_dataset(args.inFileName)

x = dsIn.x.values
y = dsIn.y.values
dx = int((x[1]-x[0])/1000.)
Lx = int((x[-1] - x[0])/1000.)
Ly = int((y[-1] - y[0])/1000.)

inMeshName = '{}x{}km_{}km_Antarctic_stereo'.format(Lx, Ly, dx)

projection = get_antarctic_stereographic_projection()

inDescriptor = ProjectionGridDescriptor.create(projection, x, y, inMeshName)

outRes = args.resolution*1e3

nxOut = int((x[-1] - x[0])/outRes + 0.5) + 1
nyOut = int((y[-1] - y[0])/outRes + 0.5) + 1

xOut = x[0] + outRes*numpy.arange(nxOut)
yOut = y[0] + outRes*numpy.arange(nyOut)


outMeshName = '{}x{}km_{}km_Antarctic_stereo'.format(Lx, Ly, args.resolution)

outDescriptor = ProjectionGridDescriptor.create(projection, xOut, yOut,
                                                outMeshName)

mappingFileName = 'map_{}_to_{}.nc'.format(inMeshName, outMeshName)

remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

remapper.build_mapping_file(method='bilinear')

dsOut = remapper.remap(dsIn, renormalizationThreshold=0.01)
dsOut.to_netcdf(args.outFileName)
