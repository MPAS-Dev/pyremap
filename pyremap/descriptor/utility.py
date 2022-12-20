# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/pyremap/main/LICENSE

import numpy


def interp_extrap_corner(inField):
    """Interpolate/extrapolate a 1D field from grid centers to grid corners"""

    outField = numpy.zeros(len(inField) + 1)
    outField[1:-1] = 0.5 * (inField[0:-1] + inField[1:])
    # extrapolate the ends
    outField[0] = 1.5 * inField[0] - 0.5 * inField[1]
    outField[-1] = 1.5 * inField[-1] - 0.5 * inField[-2]
    return outField


def interp_extrap_corners_2d(inField):
    """Interpolate/extrapolate a 1D field from grid centers to grid corners"""

    temp = numpy.zeros((inField.shape[0], inField.shape[1] + 1))
    temp[:, 1:-1] = 0.5 * (inField[:, 0:-1] + inField[:, 1:])
    # extrapolate the ends
    temp[:, 0] = 1.5 * inField[:, 0] - 0.5 * inField[:, 1]
    temp[:, -1] = 1.5 * inField[:, -1] - 0.5 * inField[:, -2]

    outField = numpy.zeros((inField.shape[0] + 1, inField.shape[1] + 1))
    outField[1:-1, :] = 0.5 * (temp[0:-1, :] + temp[1:, :])
    # extrapolate the ends
    outField[0, :] = 1.5 * temp[0, :] - 0.5 * temp[1, :]
    outField[-1, :] = 1.5 * temp[-1, :] - 0.5 * temp[-2, :]

    return outField


def create_scrip(outFile, grid_size, grid_corners, grid_rank, units,
                 meshName):
    """
    Given a SCRIP files, creates common variables and writes common values used
    in various types of SCRIP files.

    Parameters
    ----------
    outFile : file pointer
        A SCRIP file opened in write mode

    grid_size : int
        The number of elements in the grid or mesh

    grid_corners : int
        The number of corners in the grid or mesh

    grid_rank : int
        The dimensionality of the grid (1 for mesh, 2 for grid)

    units : {'degrees', 'radians'}
        The units for latitude and longitude

    meshName : str
        The name of the mesh
    """
    # Write to output file
    # Dimensions
    outFile.createDimension("grid_size", grid_size)
    outFile.createDimension("grid_corners", grid_corners)
    outFile.createDimension("grid_rank", grid_rank)

    # Variables
    grid_center_lat = outFile.createVariable('grid_center_lat', 'f8',
                                             ('grid_size',))
    grid_center_lat.units = units
    grid_center_lon = outFile.createVariable('grid_center_lon', 'f8',
                                             ('grid_size',))
    grid_center_lon.units = units
    grid_corner_lat = outFile.createVariable('grid_corner_lat', 'f8',
                                             ('grid_size', 'grid_corners'))
    grid_corner_lat.units = units
    grid_corner_lon = outFile.createVariable('grid_corner_lon', 'f8',
                                             ('grid_size', 'grid_corners'))
    grid_corner_lon.units = units
    grid_imask = outFile.createVariable('grid_imask', 'i4', ('grid_size',))
    grid_imask.units = 'unitless'
    outFile.createVariable('grid_dims', 'i4', ('grid_rank',))

    outFile.meshName = meshName


def unwrap_corners(inField):
    """Turn a 2D array of corners into an array of rectangular mesh elements"""
    outField = numpy.zeros(((inField.shape[0] - 1) * (inField.shape[1] - 1),
                            4))
    # corners are counterclockwise
    outField[:, 0] = inField[0:-1, 0:-1].flat
    outField[:, 1] = inField[0:-1, 1:].flat
    outField[:, 2] = inField[1:, 1:].flat
    outField[:, 3] = inField[1:, 0:-1].flat

    return outField


def round_res(res):
    """Round the resoltuion to a reasonable number for grid names"""
    rounded = numpy.round(res*1000.)/1000.
    return '{}'.format(rounded)
