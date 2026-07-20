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
Unit tests for the commands pyremap builds to make mapping files.
"""

import os
import shutil
import tempfile
from unittest import mock

import numpy as np
import pytest

from pyremap import LatLonGridDescriptor, Remapper
from tests import TestCase


def _get_descriptors():
    src_descriptor = LatLonGridDescriptor.create(
        lat_corner=np.linspace(-90.0, 90.0, 46),
        lon_corner=np.linspace(-180.0, 180.0, 91),
        units='degrees',
    )
    dst_descriptor = LatLonGridDescriptor.create(
        lat_corner=np.linspace(-90.0, 90.0, 91),
        lon_corner=np.linspace(-180.0, 180.0, 181),
        units='degrees',
    )
    return src_descriptor, dst_descriptor


def _build_map_calls(ntasks, map_tool='moab', method='bilinear'):
    """
    Build a map with ``check_call`` mocked out, returning the list of argument
    lists that would have been passed to ``check_call``
    """
    src_descriptor, dst_descriptor = _get_descriptors()
    remapper = Remapper(
        ntasks=ntasks,
        map_tool=map_tool,
        method=method,
        src_descriptor=src_descriptor,
        dst_descriptor=dst_descriptor,
        use_tmp=False,
    )
    with mock.patch(
        'pyremap.remapper.build_map.check_call'
    ) as mock_check_call:
        remapper.build_map()

    return [call.args[0] for call in mock_check_call.call_args_list]


class TestBuildMapArgs(TestCase):
    """
    Test the argument lists passed to the mapping tools, with ``check_call``
    mocked out so no ESMF or MOAB binaries are needed.
    """

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.test_dir)

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.test_dir)

    def test_moab_serial_does_not_partition(self):
        """
        ``mbpart`` cannot make a single partition, so it must not be called
        when ``ntasks=1``
        """
        calls = _build_map_calls(ntasks=1)

        executables = [os.path.basename(args[0]) for args in calls]
        assert 'mbpart' not in executables

        # the SCRIP files are still converted to HDF5, which mbtempest needs
        assert executables == ['mbconvert', 'mbconvert', 'mbtempest']

        mbtempest_args = calls[-1]
        # no MPI launcher and the unpartitioned .h5m meshes are used
        assert mbtempest_args[0].endswith('mbtempest')
        loaded = [
            mbtempest_args[index + 1]
            for index, arg in enumerate(mbtempest_args)
            if arg == '--load'
        ]
        assert loaded == ['src_mesh.h5m', 'dst_mesh.h5m']

    def test_moab_parallel_partitions(self):
        """
        With more than one task, the meshes are partitioned with ``mbpart``
        and mbtempest is launched with an MPI launcher
        """
        ntasks = 2
        calls = _build_map_calls(ntasks=ntasks)

        executables = [os.path.basename(args[0]) for args in calls]
        assert executables == [
            'mbconvert',
            'mbpart',
            'mbconvert',
            'mbpart',
            'mpirun',
        ]

        for args in calls:
            if os.path.basename(args[0]) == 'mbpart':
                assert args[1] == f'{ntasks}'

        mbtempest_args = calls[-1]
        assert mbtempest_args[:3] == ['mpirun', '-np', f'{ntasks}']
        loaded = [
            mbtempest_args[index + 1]
            for index, arg in enumerate(mbtempest_args)
            if arg == '--load'
        ]
        assert loaded == [
            f'src_mesh.p{ntasks}.h5m',
            f'dst_mesh.p{ntasks}.h5m',
        ]

    def test_esmf_serial(self):
        """
        The ESMF branch never partitions, whatever the task count
        """
        calls = _build_map_calls(ntasks=1, map_tool='esmf')

        assert len(calls) == 1
        assert os.path.basename(calls[0][0]) == 'ESMF_RegridWeightGen'


@pytest.mark.parametrize('ntasks', [1, 2])
def test_moab_serial_args_match_parallel(ntasks):
    """
    Whatever the task count, mbtempest gets the same weight-generation
    options
    """
    cwd = os.getcwd()
    test_dir = tempfile.mkdtemp()
    os.chdir(test_dir)
    try:
        calls = _build_map_calls(ntasks=ntasks)
    finally:
        os.chdir(cwd)
        shutil.rmtree(test_dir)

    mbtempest_args = calls[-1]
    start = mbtempest_args.index('--file')
    assert mbtempest_args[start:] == [
        '--file',
        mbtempest_args[start + 1],
        '--weights',
        '--gnomonic',
        '--boxeps',
        '1e-9',
        '--method',
        'fv',
        '--method',
        'fv',
        '--order',
        '1',
        '--order',
        '1',
        '--fvmethod',
        'bilin',
    ]
