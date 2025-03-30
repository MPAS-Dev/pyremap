#!/usr/bin/env python

"""
Creates a mapping file and remaps the variables from an input file on an
Antarctic grid to another grid with the same extent but a different resolution.
The mapping file can be used with ncremap (NCO) to remap MPAS files between
these same gridss.

Usage: Copy this script into the main MPAS-Analysis directory (up one level).
"""

import argparse

import numpy as np
import xarray as xr

from pyremap import Remapper, ProjectionGridDescriptor
from pyremap.polar import get_antarctic_stereographic_projection


parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', dest='in_filename', required=True, type=str,
                    help='Input file name')
parser.add_argument('-o', dest='out_filename', required=True, type=str,
                    help='Output file name')
parser.add_argument('-r', dest='resolution', required=True, type=float,
                    help='Output resolution')
parser.add_argument('-m', dest='method', required=False, default='bilinear',
                    help='Method: {"bilinear", "neareststod", "conserve"}')
parser.add_argument('-t', dest='mpi_tasks', required=False, type=int,
                    default=1,
                    help='Number of MPI tasks (default = 1)')
args = parser.parse_args()


if args.method not in ['bilinear', 'neareststod', 'conserve']:
    raise ValueError(f'Unexpected method {args.method}')

ds_in = xr.open_dataset(args.in_filename)

x = ds_in.x.values
y = ds_in.y.values
dx = int((x[1] - x[0]) / 1000.0)
lx = int((x[-1] - x[0]) / 1000.0)
ly = int((y[-1] - y[0]) / 1000.0)

src_mesh_name = f'{lx}x{ly}km_{dx}km_Antarctic_stereo'

projection = get_antarctic_stereographic_projection()

remapper = Remapper(ntasks=args.mpi_tasks, method=args.method)
remapper.src_descriptor = ProjectionGridDescriptor.create(
    projection, x, y, src_mesh_name)

out_res = args.resolution * 1e3

nx_out = int((x[-1] - x[0]) / out_res + 0.5) + 1
ny_out = int((y[-1] - y[0]) / out_res + 0.5) + 1

x_out = x[0] + out_res * np.arange(nx_out)
y_out = y[0] + out_res * np.arange(ny_out)

dst_mesh_name = f'{lx}x{ly}km_{args.resolution}km_Antarctic_stereo'

remapper.dst_descriptor = ProjectionGridDescriptor.create(
    projection, x_out, y_out, dst_mesh_name)

remapper.build_map()

ds_out = remapper.remap_numpy(ds_in, renormalization_threshold=0.01)
ds_out.to_netcdf(args.out_filename)
