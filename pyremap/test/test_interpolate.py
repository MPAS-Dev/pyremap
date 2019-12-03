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
    PointCollectionDescriptor, write_netcdf

# noinspection PyUnresolvedReferences
from pyremap.test import TestCase, loaddatadir


# noinspection PyMethodMayBeStatic
@pytest.mark.usefixtures('loaddatadir')
class TestInterp(TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)
        #pass

    def get_mpas_descriptor(self):
        mpas_mesh_filename = str(self.datadir.join('mpasMesh.nc'))
        time_series_filename = \
            str(self.datadir.join('timeSeries.0002-01-01.nc'))
        outfilename = '{}/unmapped_mpas.nc'.format(self.test_dir)
        # add fill values
        write_netcdf(xarray.open_dataset(time_series_filename),outfilename,
                     always_fill=True)
        descriptor = MpasMeshDescriptor(mpas_mesh_filename, meshname='oQU240')

        return descriptor, mpas_mesh_filename, outfilename

    def get_lon_lat_file_descriptor(self):
        lon_lat_grid_filename = \
            str(self.datadir.join('SST_annual_1870-1900.nc'))
        descriptor = LonLatGridDescriptor(filename=lon_lat_grid_filename,
                                          lonvarname='lon',
                                          latvarname='lat',
                                          units='degrees')

        return descriptor, lon_lat_grid_filename

    def get_lon_lat_2d_descriptor(self):
        infilename = str(self.datadir.join('stereographic_test.nc'))

        meshname = '100km_Antarctic_stereo'
        descriptor = LonLat2DGridDescriptor(filename=infilename,
                                            meshname=meshname,
                                            lonvarname='lon',
                                            latvarname='lat',
                                            units='degrees')

        return descriptor, infilename

    def get_constant_spacing_lon_lat_descriptor(self, res=2.):
        descriptor = LonLatGridDescriptor.create_constant_spacing(
            dlon=res, dlat=res, lon_min=-180., lon_max=180., lat_min=-90.,
            lat_max=90.)
        return descriptor

    def get_stereographic_file_descriptor(self):

        # projection for BEDMAP2 and other common Antarctic data sets
        projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 '
                                 '+lon_0=0.0  +k_0=1.0 +x_0=0.0 +y_0=0.0 '
                                 '+ellps=WGS84')

        infilename = str(self.datadir.join('stereographic_test.nc'))
        meshname = '100km_Antarctic_stereo'
        descriptor = ProjectionGridDescriptor(projection=projection,
                                              filename=infilename,
                                              meshname=meshname)
        return descriptor, infilename


    def get_stereographic_array_descriptor(self, res=100e3):

        # projection for BEDMAP2 and other common Antarctic data sets
        projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 '
                                 '+lon_0=0.0  +k_0=1.0 +x_0=0.0 +y_0=0.0 '
                                 '+ellps=WGS84')

        # a square 61x61 cell map with 100 km resolution and
        xmax = 3000e3
        nx = 2 * int(xmax / res) + 1
        x = numpy.linspace(-xmax, xmax, nx)
        meshname = '{}km_Antarctic_stereo'.format(int(res * 1e-3))
        descriptor = ProjectionGridDescriptor(projection=projection, x=x, y=x,
                                              meshname=meshname, units='meters')
        return descriptor

    def get_points_descriptor(self, npoints=1000):

        numpy.random.seed(seed=0)
        lon = 360.*numpy.random.rand(npoints) - 180.
        lat = 180.*numpy.random.rand(npoints) - 90.
        meshname = 'random_points'
        descriptor = PointCollectionDescriptor(
            lats=lat, lons=lon, collection_name=meshname, units='degrees',
            out_dimension='npoints')
        return descriptor

    def get_file_names(self, suffix):
        mapping_filename = '{}/weights_{}.nc'.format(self.test_dir, suffix)
        outfilename = '{}/remapped_{}.nc'.format(self.test_dir, suffix)
        ref_filename = '{}/ref_{}.nc'.format(self.datadir, suffix)
        return mapping_filename, outfilename, ref_filename

    def build_remapper(self, src_descrip, dst_descrip, mapping_filename,
                       method='bilinear'):

        remapper = Remapper(src_descrip, dst_descrip, mapping_filename)

        remapper.build_mapping_file(method=method)

        assert os.path.exists(remapper.mapping_filename)

        return remapper

    def check_remap(self, infilename, outfilename, ref_filename, remapper,
                    remap_file=True, renorm=0.01):

        if os.path.exists(ref_filename):
            ds_ref = self.drop_extras(xarray.open_dataset(ref_filename))
        else:
            ds = xarray.open_dataset(infilename)
            ds_ref = remapper.remap(
                ds=ds, renorm_thresh=renorm)
            write_netcdf(ds_ref, os.path.basename(ref_filename))

        if remap_file:
            # first, test interpolation with ncremap
            remapper.remap_file(infilename=infilename,
                                outfilename=outfilename,
                                renormalize=renorm)

            assert os.path.exists(outfilename)
            ds_remapped = self.drop_extras(xarray.open_dataset(outfilename))
            print(ds_remapped.data_vars, ds_ref.data_vars)
            print(ds_remapped.coords, ds_ref.coords)
            self.assertDatasetApproxEqual(ds_remapped, ds_ref)

        # now, try in-memory remapping
        ds = xarray.open_dataset(infilename)
        ds_remapped = remapper.remap(
            ds=ds, renorm_thresh=renorm)
        ds_remapped = self.drop_extras(ds_remapped)
        self.assertDatasetApproxEqual(ds_remapped, ds_ref)

    def drop_extras(self, ds, extra_vars=None):
        if extra_vars is None:
            extra_vars = ['lat_bnds', 'lon_bnds', 'gw', 'area', 'x_bnds',
                          'y_bnds', 'lat_vertices', 'lon_vertices']
        varnames = [var for var in extra_vars if var in ds]
        ds = ds.drop_vars(varnames)
        return ds

    def make_stereographic_dataset(self, src_descrip):
        Lat = src_descrip.coords['lat']['data']

        # now, let's make a more complicated field with more dimensions to
        # check
        infield = numpy.reshape(Lat, (1, Lat.shape[0], Lat.shape[1], 1))
        infield = infield.repeat(3, axis=0)
        infield = infield.repeat(2, axis=3)

        all_ones = numpy.ones(Lat.shape)

        datasetDict = {'dims': ('dim0', 'y', 'x', 'dim3'),
                       'coords': src_descrip.coords,
                       'data_vars': {'complicated':
                                     {'dims': ('dim0', 'y', 'x', 'dim3'),
                                      'data': infield},
                                     'all_ones':
                                         {'dims': ('y', 'x'),
                                          'data': all_ones}}}

        ds = xarray.Dataset.from_dict(datasetDict)
        infilename = '{}/unmapped_stereographic.nc'.format(self.test_dir)
        write_netcdf(ds, infilename)
        return infilename

    def make_lon_lat_dataset(self, src_descrip):

        Lon, Lat = numpy.meshgrid(src_descrip.coords['lon']['data'],
                                  src_descrip.coords['lat']['data'])
        infield = Lat**2

        datasetDict = {'dims': ('lat', 'lon'),
                       'coords': src_descrip.coords,
                       'data_vars': {'lat_squared':
                                     {'dims': ('lat', 'lon'),
                                      'data': infield}}}

        ds = xarray.Dataset.from_dict(datasetDict)
        infilename = '{}/unmapped_lon_lat.nc'.format(self.test_dir)
        write_netcdf(ds, infilename)
        return infilename

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

        infilename = self.make_stereographic_dataset(src_descrip)

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)

        self.check_remap(infilename, outfilename, ref_filename,
                         remapper, remap_file=True)

    def test_lat_lon_2d_to_lat_lon(self):
        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='lon_lat_2d_to_lon_lat')

        src_descrip, infilename = self.get_lon_lat_2d_descriptor()
        dst_descrip = self.get_constant_spacing_lon_lat_descriptor()

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)

        self.check_remap(infilename, outfilename, ref_filename,
                         remapper, remap_file=True)

    def test_stereographic_to_stereographic(self):
        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='stereographic_to_stereographic')

        src_descrip, infilename = self.get_stereographic_file_descriptor()
        dst_descrip = self.get_stereographic_array_descriptor(res=200e3)

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)
        self.check_remap(infilename, outfilename, ref_filename,
                         remapper, remap_file=True)

    def test_mpas_to_stereographic_nearest(self):
        """
        test horizontal interpolation from an MPAS mesh to a destination
        stereographic grid.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='mpas_to_stereographic_nearest')

        src_descrip, mpas_mesh_filename, time_series_filename = \
            self.get_mpas_descriptor()
        dst_descrip, _ = self.get_stereographic_file_descriptor()

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename,
                                       method='neareststod')

        self.check_remap(time_series_filename, outfilename, ref_filename,
                         remapper, remap_file=True)

    # noinspection PyTypeChecker
    def test_lon_lat_to_stereographic_conserve(self):
        """
        test horizontal interpolation from an MPAS mesh to a destination
        stereographic grid.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='lon_lat_to_stereographic_conserve')

        src_descrip = self.get_constant_spacing_lon_lat_descriptor(res=10.)
        dst_descrip = self.get_stereographic_array_descriptor(res=500e3)

        infilename = self.make_lon_lat_dataset(src_descrip)

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename, method='conserve')

        self.check_remap(infilename, outfilename, ref_filename,
                         remapper, remap_file=True, renorm=None)

    # noinspection PyTypeChecker
    def test_stereographic_to_lon_lat_conserve(self):
        """
        test horizontal interpolation from an MPAS mesh to a destination
        stereographic grid.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='stereographic_to_lon_lat_conserve')

        src_descrip = self.get_stereographic_array_descriptor(res=500e3)
        dst_descrip = self.get_constant_spacing_lon_lat_descriptor(res=10.)

        infilename = self.make_stereographic_dataset(src_descrip)

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename, method='conserve')

        self.check_remap(infilename, outfilename, ref_filename,
                         remapper, remap_file=True, renorm=None)

    def test_lon_lat_to_points(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        stereographic grid.
        """

        mapping_filename, outfilename, ref_filename = \
            self.get_file_names(suffix='lon_lat_to_points')

        src_descrip, lon_lat_grid_filename = self.get_lon_lat_file_descriptor()
        dst_descrip = self.get_points_descriptor()

        remapper = self.build_remapper(src_descrip, dst_descrip,
                                       mapping_filename)

        self.check_remap(lon_lat_grid_filename, outfilename, ref_filename,
                         remapper, remap_file=True)
