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
Unit test infrastructure for horizontal interpolation.
"""

import pytest
import shutil
import os
import tempfile
import numpy
import xarray
import pyproj

from pyremap import Remapper, MpasMeshDescriptor, \
    LonLatGridDescriptor, ProjectionGridDescriptor, LonLat2DGridDescriptor, \
    PointCollectionDescriptor, MeshDescriptor, write_netcdf

from pyremap.test import TestCase, loaddatadir


# noinspection PyMethodMayBeStatic
@pytest.mark.usefixtures('loaddatadir')
class TestInterp(TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        self.renormalizationThreshold = 0.01

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)
        pass

    def get_mpas_descriptor(self):
        mpas_mesh_filename = str(self.datadir.join('mpasMesh.nc'))
        time_series_filename = \
            str(self.datadir.join('timeSeries.0002-01-01.nc'))
        descriptor = MpasMeshDescriptor(mpas_mesh_filename, meshname='oQU240')

        return descriptor, mpas_mesh_filename, time_series_filename

    def get_lon_lat_file_descriptor(self):
        lon_lat_grid_filename = \
            str(self.datadir.join('SST_annual_1870-1900.nc'))
        descriptor = LonLatGridDescriptor(filename=lon_lat_grid_filename,
                                          lonvarname='lon',
                                          latvarname='lat',
                                          units='degrees')

        return descriptor, lon_lat_grid_filename

    def get_constant_spacing_lon_lat_descriptor(self):
        descriptor = LonLatGridDescriptor.create_constant_spacing(
            dlon=2., dlat=2., lon_min=-180., lon_max=180., lat_min=-90.,
            lat_max=90.)
        return descriptor

    def get_stereographic_array_descriptor(self):

        # projection for BEDMAP2 and other common Antarctic data sets
        projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 '
                                 '+lon_0=0.0  +k_0=1.0 +x_0=0.0 +y_0=0.0 '
                                 '+ellps=WGS84')

        # a square 61x61 cell map with 100 km resolution and
        xmax = 3000e3
        res = 100e3
        nx = 2 * int(xmax / res) + 1
        x = numpy.linspace(-xmax, xmax, nx)
        meshname = '{}km_Antarctic_stereo'.format(int(res * 1e-3))
        descriptor = ProjectionGridDescriptor(projection=projection, x=x, y=x,
                                              meshname=meshname, units='meters')
        return descriptor

    def get_file_names(self, suffix):
        mapping_filename = '{}/weights_{}.nc'.format(self.test_dir, suffix)
        outfilename = '{}/remapped_{}.nc'.format(self.test_dir, suffix)
        ref_filename = '{}/ref_{}.nc'.format(self.datadir, suffix)
        return mapping_filename, outfilename, ref_filename

    @staticmethod
    def build_remapper(src_descrip, dst_descrip, mapping_filename):

        remapper = Remapper(src_descrip, dst_descrip, mapping_filename)

        remapper.build_mapping_file(method='bilinear')

        assert os.path.exists(remapper.mapping_filename)

        return remapper

    def check_remap(self, infilename, outfilename, ref_filename, remapper,
                    remap_file=True):

        ds_ref = xarray.open_dataset(ref_filename)
        if remap_file:
            # first, test interpolation with ncremap
            remapper.remap_file(infilename=infilename,
                                outfilename=outfilename)

            assert os.path.exists(outfilename)
            ds_remapped = xarray.open_dataset(outfilename)
            # drop some extra variables added by ncremap that aren't in the
            # reference data set
            varnames = [var for var in ['lat_bnds', 'lon_bnds', 'gw', 'area']
                        if var in ds_remapped]
            ds_remapped = ds_remapped.drop_vars(varnames)
            self.assertDatasetApproxEqual(ds_remapped, ds_ref)

        # now, try in-memory remapping
        ds = xarray.open_dataset(infilename)
        ds_remapped = remapper.remap(ds, self.renormalizationThreshold)
        self.assertDatasetApproxEqual(ds_remapped, ds_ref)

    def test_mpas_to_lon_lat_file(self):
        """
        test horizontal interpolation from an MPAS mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='mpas_to_lon_lat_file')

        src_descrip, mpas_mesh_filename, time_series_filename = \
            self.get_mpas_descriptor()
        dst_descrip, lon_lat_grid_filename = \
            self.get_lon_lat_file_descriptor()

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)
        self.check_remap(time_series_filename, outfilename, ref_filename,
                         remapper, remap_file=True)

    def test_mpas_to_lon_lat_array(self):
        """
        test horizontal interpolation from an MPAS mesh to a destination
        lat/lon grid determined from config options 'lat' and 'lon'.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='mpas_to_lon_lat_array')

        src_descrip, mpas_mesh_filename, time_series_filename = \
            self.get_mpas_descriptor()
        dst_descrip = self.get_constant_spacing_lon_lat_descriptor()

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)
        self.check_remap(time_series_filename, outfilename, ref_filename,
                         remapper, remap_file=True)

    def test_lon_lat_file_to_lon_lat_array(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        lat/lon grid determined from config options 'lat' and 'lon'.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='lon_lat_file_to_lon_lat_array')

        src_descrip, lon_lat_grid_filename = \
            self.get_lon_lat_file_descriptor()
        dst_descrip = self.get_constant_spacing_lon_lat_descriptor()

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)
        self.check_remap(lon_lat_grid_filename, outfilename, ref_filename,
                         remapper, remap_file=True)

    def test_mpas_to_stereographic_array(self):
        """
        test horizontal interpolation from an MPAS mesh to a destination
        stereographic grid.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='mpas_to_stereographic_array')

        src_descrip, mpas_mesh_filename, time_series_filename = \
            self.get_mpas_descriptor()
        dst_descrip = self.get_stereographic_array_descriptor()

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)

        self.check_remap(time_series_filename, outfilename, ref_filename,
                         remapper, remap_file=True)

    def test_lon_lat_file_to_stereographic_array(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        stereographic grid.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='lon_lat_file_to_stereographic_array')

        src_descrip, lon_lat_grid_filename = \
            self.get_lon_lat_file_descriptor()
        dst_descrip = self.get_stereographic_array_descriptor()

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)

        self.check_remap(lon_lat_grid_filename, outfilename, ref_filename,
                         remapper, remap_file=True)

    def test_stereographic_array_to_lon_lat_array(self):
        """
        test horizontal interpolation from a stereographic grid to a
        destination lat/lon grid determined from config options 'lat' and
        'lon'.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='stereographic_array_to_lon_lat_array')

        src_descrip = self.get_stereographic_array_descriptor()
        dst_descrip = self.get_constant_spacing_lon_lat_descriptor()

        Lat = src_descrip.coords['lat']['data']

        # now, let's make a more complicated field with more dimensions to
        # check
        inField = numpy.reshape(Lat, (1, Lat.shape[0], Lat.shape[1], 1))
        inField = inField.repeat(3, axis=0)
        inField = inField.repeat(2, axis=3)

        datasetDict = {'dims': ('dim0', 'x', 'y', 'dim3'),
                       'coords': src_descrip.coords,
                       'data_vars': {'complicated':
                                     {'dims': ('dim0', 'x', 'y', 'dim3'),
                                      'data': inField}}}

        ds = xarray.Dataset.from_dict(datasetDict)
        infilename = '{}/unmapped_stereographic_array_to_lon_lat_' \
                      'array.nc'.format(self.test_dir)
        write_netcdf(ds, infilename)

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                             mapping_filename)

        self.check_remap(infilename, outfilename, ref_filename,
                         remapper, remap_file=True)
