#!/usr/bin/env python

import xarray
import numpy as np
import os
import subprocess
from importlib.resources import path


def check_remap(outFileName, refFileName):

    dsRef = xarray.open_dataset(refFileName)
    assert os.path.exists(outFileName)
    dsRemapped = xarray.open_dataset(outFileName)
    # drop some extra vairables added by ncremap that aren't in the
    # reference data set
    dsRemapped = dsRemapped.drop_vars(['lat_bnds', 'lon_bnds', 'gw',
                                       'area'])
    assertDatasetApproxEqual(dsRemapped, dsRef)


def assertArrayApproxEqual(a1, a2, rtol=1e-5, atol=1e-8):
    close = np.isclose(a1, a2, rtol=rtol, atol=atol)
    oneIsNaN = np.logical_or(np.isnan(a1), np.isnan(a2))
    assert np.all(np.logical_or(close, oneIsNaN))


def assertDatasetApproxEqual(ds1, ds2, rtol=1e-5, atol=1e-8):
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
            assertArrayApproxEqual(ds1[var].values, ds2[var].values,
                                   rtol=rtol, atol=atol)
    else:
        assertArrayApproxEqual(ds1.values, ds2.values, rtol=rtol, atol=atol)


def test_ncremap():
    test_dir = 'test_ncremap'
    ref_package = 'pyremap.test.ncremap_files'

    try:
        os.makedirs(test_dir)
    except OSError:
        pass

    with path(ref_package, 'timeSeries.0002-01-01.nc') as filename:
        inFileName = str(filename)
    with path(ref_package, 'weights_mpas_to_latlon_file.nc') as filename:
        weightFileName = str(filename)
    with path(ref_package, 'ref_mpas_to_latlon_file.nc') as filename:
        refFileName = str(filename)
    outFileName = os.path.join(test_dir, 'remapped_mpas_to_latlon_file.nc')

    command = 'ncremap -m {} --vrb=1 ' \
              '-R "--rgr lat_nm_out=lat --rgr lon_nm_out=lon" ' \
              '-P mpas -C {} {}'.format(weightFileName, inFileName,
                                        outFileName)

    print('running: {}'.format(command))
    subprocess.check_call(command, shell=True)

    check_remap(outFileName, refFileName)
