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
Unit test infrastructure, adapted from approach of xarray.
"""

import os
import unittest
import warnings
from contextlib import contextmanager
from shutil import copytree, ignore_patterns

import numpy as np
import xarray
from pytest import fixture


# Adapted from
# http://stackoverflow.com/questions/29627341/pytest-where-to-store-expected-data
@fixture
def loaddatadir(request, tmpdir):
    """
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    """
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)
    print('tmpdir: ', tmpdir)
    print('test_dir: ', test_dir)

    if os.path.isdir(test_dir):
        copytree(test_dir, str(tmpdir), dirs_exist_ok=True,
                 ignore=ignore_patterns('__pycache__'))

    request.cls.datadir = tmpdir


class TestCase(unittest.TestCase):
    def assertEqual(self, a1, a2):
        assert a1 == a2 or (a1 != a1 and a2 != a2)

    def assertLessThan(self, a1, a2):
        assert a1 <= a2

    def assertGreaterThan(self, a1, a2):
        assert a1 >= a2

    def assertArrayEqual(self, a1, a2):
        np.testing.assert_array_equal(a1, a2)

    def assertApproxEqual(self, a1, a2, rtol=1e-5, atol=1e-8):
        assert np.isclose(a1, a2, rtol=rtol, atol=atol)

    def assertArrayApproxEqual(self, a1, a2, rtol=1e-5, atol=1e-8):
        close = np.isclose(a1, a2, rtol=rtol, atol=atol)
        oneIsNaN = np.logical_or(np.isnan(a1), np.isnan(a2))
        assert np.all(np.logical_or(close, oneIsNaN))

    def assertDatasetApproxEqual(self, ds1, ds2, rtol=1e-5, atol=1e-8):
        assert ((isinstance(ds1, xarray.Dataset) and
                 isinstance(ds2, xarray.Dataset)) or
                (isinstance(ds1, xarray.DataArray) and
                 isinstance(ds2, xarray.DataArray)))

        if isinstance(ds1, xarray.Dataset):
            assert set(ds1.data_vars.keys()) == set(ds2.data_vars.keys())
            for var in ds1.data_vars.keys():
                assert var in ds2.data_vars.keys()
                if ds1[var].values.dtype.char == 'S':
                    continue
                self.assertArrayApproxEqual(ds1[var].values, ds2[var].values,
                                            rtol=rtol, atol=atol)
        else:
            self.assertArrayApproxEqual(ds1.values, ds2.values,
                                        rtol=rtol, atol=atol)

    @contextmanager
    def assertWarns(self, message):
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', message)
            yield
            assert len(w) > 0
            assert all(message in str(wi.message) for wi in w)
