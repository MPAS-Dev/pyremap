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
Unit test infrastructure for horizontal interpolation.
"""

import os
import shutil
import tempfile
from typing import Optional

import numpy as np
import pyproj
import pytest
import xarray as xr

from pyremap import (
    LatLon2DGridDescriptor,
    LatLonGridDescriptor,
    MpasCellMeshDescriptor,
    MpasEdgeMeshDescriptor,
    MpasVertexMeshDescriptor,
    PointCollectionDescriptor,
    ProjectionGridDescriptor,
    Remapper,
)
from pyremap.test import TestCase, loaddatadir  # noqa: F401


@pytest.mark.usefixtures('loaddatadir')
class TestInterp(TestCase):
    datadir: Optional[str]  # Added type hint for datadir

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        self.renormalization_threshold = 0.01
        self.cwd = os.getcwd()
        os.chdir(self.test_dir)

    def tearDown(self):
        # Remove the directory after the test
        os.chdir(self.cwd)
        shutil.rmtree(self.test_dir)

    def get_mpas_cell_descriptor(self):
        assert self.datadir is not None, 'datadir must not be None'
        mpas_mesh_filename = str(self.datadir.join('mpasMesh.nc'))
        time_series_filename = str(
            self.datadir.join('timeSeries.0002-01-01.nc')
        )
        descriptor = MpasCellMeshDescriptor(
            mpas_mesh_filename, mesh_name='oQU240'
        )

        return descriptor, mpas_mesh_filename, time_series_filename

    def get_mpas_edge_descriptor(self):
        assert self.datadir is not None, 'datadir must not be None'
        mpas_mesh_filename = str(self.datadir.join('mpasMesh.nc'))
        area_filename = str(self.datadir.join('mpasAreaEdge.nc'))
        descriptor = MpasEdgeMeshDescriptor(
            mpas_mesh_filename, mesh_name='oQU240'
        )

        return descriptor, mpas_mesh_filename, area_filename

    def get_mpas_vertex_descriptor(self):
        assert self.datadir is not None, 'datadir must not be None'
        mpas_mesh_filename = str(self.datadir.join('mpasMesh.nc'))
        area_filename = str(self.datadir.join('mpasAreaVertex.nc'))
        descriptor = MpasVertexMeshDescriptor(
            mpas_mesh_filename, mesh_name='oQU240'
        )

        return descriptor, mpas_mesh_filename, area_filename

    def get_latlon_file_descriptor(self):
        assert self.datadir is not None, 'datadir must not be None'
        latlon_grid_filename = str(
            self.datadir.join('SST_annual_1870-1900.nc')
        )
        descriptor = LatLonGridDescriptor.read(
            latlon_grid_filename, lat_var_name='lat', lon_var_name='lon'
        )

        return descriptor, latlon_grid_filename

    def get_latlon_array_descriptor(self):
        lat = np.linspace(-90.0, 90.0, 91)
        lon = np.linspace(-180.0, 180.0, 181)

        descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')
        return descriptor

    def get_latlon2d_file_descriptor(self):
        assert self.datadir is not None, 'datadir must not be None'
        latlon_grid_filename = str(
            self.datadir.join('SST_annual_1870-1900.nc')
        )
        ds = xr.open_dataset(latlon_grid_filename)
        lon2d, lat2d = np.meshgrid(ds.lon.values, ds.lat.values)
        ds['lat2d'] = (('lat', 'lon'), lat2d)
        ds['lon2d'] = (('lat', 'lon'), lon2d)
        ds.lat2d.attrs['units'] = ds.lat.attrs['units']
        ds.lon2d.attrs['units'] = ds.lon.attrs['units']
        descriptor = LatLon2DGridDescriptor.read(
            ds=ds, lat_var_name='lat2d', lon_var_name='lon2d'
        )

        return descriptor, latlon_grid_filename

    def get_point_collection_descriptor(self):
        assert self.datadir is not None, 'datadir must not be None'
        mpas_mesh_filename = str(self.datadir.join('mpasMesh.nc'))
        ds = xr.open_dataset(mpas_mesh_filename)
        lats = ds.latCell.values
        lons = ds.lonCell.values

        descriptor = PointCollectionDescriptor(
            lats=lats,
            lons=lons,
            collection_name='mpasCellCenters',
            units='radians',
        )
        return descriptor

    def get_stereographic_array_descriptor(self):
        # projection for BEDMAP2 and other common Antarctic data sets
        projection = pyproj.Proj(
            '+proj=stere +lat_ts=-71.0 +lat_0=-90 '
            '+lon_0=0.0  +k_0=1.0 +x_0=0.0 +y_0=0.0 '
            '+ellps=WGS84'
        )

        # a 61x51 cell map with 100 km resolution and
        x_max = 3000e3
        y_max = 2500e3
        res = 100e3
        nx = 2 * int(x_max / res) + 1
        ny = 2 * int(y_max / res) + 1
        x = np.linspace(-x_max, x_max, nx)
        y = np.linspace(-y_max, y_max, ny)
        mesh_name = f'{int(res * 1e-3)}km_Antarctic_stereo'
        descriptor = ProjectionGridDescriptor.create(
            projection, x, y, mesh_name
        )
        return descriptor

    def get_filenames(self, suffix):
        assert self.datadir is not None, 'datadir must not be None'
        weight_filename = f'{self.test_dir}/weights_{suffix}.nc'
        out_filename = f'{self.test_dir}/remapped_{suffix}.nc'
        ref_filename = f'{self.datadir}/ref_{suffix}.nc'
        return weight_filename, out_filename, ref_filename

    def get_scrip_filenames(self, suffix):
        assert self.datadir is not None, 'datadir must not be None'
        scrip_filename = f'{self.test_dir}/scrip_{suffix}.nc'
        ref_filename = f'{self.datadir}/ref_scrip_{suffix}.nc'
        return scrip_filename, ref_filename

    def build_remapper(self, src_descriptor, dst_descriptor, weight_filename):
        remapper = Remapper(
            ntasks=1,
            map_filename=weight_filename,
            method='bilinear',
            use_tmp=False,
        )
        remapper.src_descriptor = src_descriptor
        remapper.dst_descriptor = dst_descriptor
        remapper.build_map()

        assert os.path.exists(remapper.map_filename)

        return remapper

    def check_scrip(self, scrip_filename, ref_filename):
        ds_ref = xr.open_dataset(ref_filename)
        assert os.path.exists(scrip_filename)
        ds_scrip = xr.open_dataset(scrip_filename)
        self.assertDatasetApproxEqual(ds_scrip, ds_ref)

    def check_remap(
        self,
        in_filename,
        out_filename,
        ref_filename,
        remapper,
        remap_file=True,
    ):
        drop_vars = [
            'lat_bnds',
            'lon_bnds',
            'gw',
            'area',
            'nvertices',
            'lat_vertices',
            'lon_vertices',
        ]
        ds_remapped_file = None
        if remap_file:
            # first, test interpolation with ncremap
            remapper.ncremap(
                in_filename=in_filename,
                out_filename=out_filename,
                replace_mpas_fill=True,
            )

            assert os.path.exists(out_filename)
            ds_remapped_file = xr.open_dataset(out_filename)
            # drop some extra variables added by ncremap that aren't in the
            # reference data set
            ds_remapped_file = ds_remapped_file.drop_vars(
                [var for var in drop_vars if var in ds_remapped_file]
            )
            new_ref_filename = out_filename.replace('.nc', '_ref.nc')
            ds_remapped_file.to_netcdf(new_ref_filename)

        # now, try in-memory remapping
        ds = xr.open_dataset(in_filename)
        ds_remapped = remapper.remap_numpy(ds, self.renormalization_threshold)

        # see how we compare with the reference data set
        ds_ref = xr.open_dataset(ref_filename)
        ds_ref = ds_ref.drop_vars([var for var in drop_vars if var in ds_ref])
        if remap_file:
            self.assertDimsEqual(ds_remapped_file, ds_ref)
            self.assertDatasetApproxEqual(ds_remapped_file, ds_ref)

        self.assertDimsEqual(ds_remapped, ds_ref)
        self.assertDatasetApproxEqual(ds_remapped, ds_ref)

    def test_latlon_file_scrip(self):
        """
        test writing a SCRIP file for a lat/lon grid file
        """

        scrip_filename, ref_filename = self.get_scrip_filenames(
            suffix='latlon_file'
        )

        src_descriptor, _ = self.get_latlon_file_descriptor()
        print(f'Writing SCRIP file for lat/lon grid {scrip_filename}')
        src_descriptor.to_scrip(scrip_filename)

        self.check_scrip(scrip_filename, ref_filename)

    def test_latlon_array_scrip(self):
        """
        test writing a SCRIP file for a lat/lon grid array
        """

        scrip_filename, ref_filename = self.get_scrip_filenames(
            suffix='latlon_array'
        )

        src_descriptor = self.get_latlon_array_descriptor()
        print(f'Writing SCRIP file for lat/lon grid {scrip_filename}')
        src_descriptor.to_scrip(scrip_filename)

        self.check_scrip(scrip_filename, ref_filename)

    def test_latlon2d_scrip(self):
        """
        test writing a SCRIP file for a 2D lat/lon grid
        """

        scrip_filename, ref_filename = self.get_scrip_filenames(
            suffix='latlon2d'
        )

        src_descriptor, _ = self.get_latlon2d_file_descriptor()
        print(f'Writing SCRIP file for 2D lat/lon grid {scrip_filename}')
        src_descriptor.to_scrip(scrip_filename)

        self.check_scrip(scrip_filename, ref_filename)

    def test_mpas_cell_scrip(self):
        """
        test writing a SCRIP file for an MPAS cell mesh
        """

        scrip_filename, ref_filename = self.get_scrip_filenames(
            suffix='mpas_cell'
        )

        src_descriptor, _, _ = self.get_mpas_cell_descriptor()
        print(f'Writing SCRIP file for MPAS cell mesh {scrip_filename}')
        src_descriptor.to_scrip(scrip_filename)

        self.check_scrip(scrip_filename, ref_filename)

    def test_mpas_edge_scrip(self):
        """
        test writing a SCRIP file for an MPAS edge mesh
        """

        scrip_filename, ref_filename = self.get_scrip_filenames(
            suffix='mpas_edge'
        )

        src_descriptor, _, _ = self.get_mpas_edge_descriptor()
        print(f'Writing SCRIP file for MPAS edge mesh {scrip_filename}')
        src_descriptor.to_scrip(scrip_filename)

        self.check_scrip(scrip_filename, ref_filename)

    def test_mpas_vertex_scrip(self):
        """
        test writing a SCRIP file for an MPAS vertex mesh
        """

        scrip_filename, ref_filename = self.get_scrip_filenames(
            suffix='mpas_vertex'
        )

        src_descriptor, _, _ = self.get_mpas_vertex_descriptor()
        print(f'Writing SCRIP file for MPAS vertex mesh {scrip_filename}')
        src_descriptor.to_scrip(scrip_filename)

        self.check_scrip(scrip_filename, ref_filename)

    def test_point_collection_scrip(self):
        """
        test writing a SCRIP file for a point collection
        """

        scrip_filename, ref_filename = self.get_scrip_filenames(
            suffix='point_collection'
        )

        src_descriptor = self.get_point_collection_descriptor()
        print(f'Writing SCRIP file for point collection {scrip_filename}')
        src_descriptor.to_scrip(scrip_filename)

        self.check_scrip(scrip_filename, ref_filename)

    def test_stereographic_scrip(self):
        """
        test writing a SCRIP file for a stereographic projection grid
        """

        scrip_filename, ref_filename = self.get_scrip_filenames(
            suffix='stereographic'
        )

        src_descriptor = self.get_stereographic_array_descriptor()
        print(f'Writing SCRIP file for stereographic grid {scrip_filename}')
        src_descriptor.to_scrip(scrip_filename)

        self.check_scrip(scrip_filename, ref_filename)

    def test_mpas_cell_to_latlon(self):
        """
        test horizontal interpolation from an MPAS cell mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords
        """

        weight_filename, out_filename, ref_filename = self.get_filenames(
            suffix='mpas_cell_to_latlon'
        )

        src_descriptor, _, time_series_filename = (
            self.get_mpas_cell_descriptor()
        )
        dst_descriptor, _ = self.get_latlon_file_descriptor()

        remapper = self.build_remapper(
            src_descriptor, dst_descriptor, weight_filename
        )
        self.check_remap(
            time_series_filename,
            out_filename,
            ref_filename,
            remapper,
            remap_file=True,
        )

    def test_mpas_edge_to_latlon(self):
        """
        test horizontal interpolation from an MPAS edge mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords
        """

        weight_filename, out_filename, ref_filename = self.get_filenames(
            suffix='mpas_edge_to_latlon'
        )

        src_descriptor, _, area_filename = self.get_mpas_edge_descriptor()
        dst_descriptor, _ = self.get_latlon_file_descriptor()

        remapper = self.build_remapper(
            src_descriptor, dst_descriptor, weight_filename
        )
        self.check_remap(
            area_filename,
            out_filename,
            ref_filename,
            remapper,
            remap_file=True,
        )

    def test_mpas_vertex_to_latlon(self):
        """
        test horizontal interpolation from an MPAS vertex mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords
        """

        weight_filename, out_filename, ref_filename = self.get_filenames(
            suffix='mpas_vertex_to_latlon'
        )

        src_descriptor, _, area_filename = self.get_mpas_vertex_descriptor()
        dst_descriptor, _ = self.get_latlon_file_descriptor()

        remapper = self.build_remapper(
            src_descriptor, dst_descriptor, weight_filename
        )
        self.check_remap(
            area_filename,
            out_filename,
            ref_filename,
            remapper,
            remap_file=True,
        )

    def test_latlon_file_to_latlon_array(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        lat/lon grid determined from config options 'lat' and 'lon'.
        """

        weight_filename, out_filename, ref_filename = self.get_filenames(
            suffix='latlon_file_to_latlon_array'
        )

        src_descriptor, latlon_grid_filename = (
            self.get_latlon_file_descriptor()
        )
        dst_descriptor = self.get_latlon_array_descriptor()

        remapper = self.build_remapper(
            src_descriptor, dst_descriptor, weight_filename
        )
        self.check_remap(
            latlon_grid_filename,
            out_filename,
            ref_filename,
            remapper,
            remap_file=True,
        )

    def test_mpas_cell_to_stereographic(self):
        """
        test horizontal interpolation from an MPAS cell mesh to a destination
        stereographic grid.
        """

        weight_filename, out_filename, ref_filename = self.get_filenames(
            suffix='mpas_cell_to_stereographic'
        )

        src_descriptor, _, time_series_filename = (
            self.get_mpas_cell_descriptor()
        )
        dst_descriptor = self.get_stereographic_array_descriptor()

        remapper = self.build_remapper(
            src_descriptor, dst_descriptor, weight_filename
        )

        # ncremap doesn't support stereographic grids so just check the
        # Remapper
        self.check_remap(
            time_series_filename,
            out_filename,
            ref_filename,
            remapper,
            remap_file=True,
        )

    def test_latlon_to_stereographic(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        stereographic grid.
        """

        weight_filename, out_filename, ref_filename = self.get_filenames(
            suffix='latlon_to_stereographic'
        )

        src_descriptor, latlon_grid_filename = (
            self.get_latlon_file_descriptor()
        )
        dst_descriptor = self.get_stereographic_array_descriptor()

        remapper = self.build_remapper(
            src_descriptor, dst_descriptor, weight_filename
        )

        # ncremap doesn't support stereographic grids so just check the
        # Remapper
        self.check_remap(
            latlon_grid_filename,
            out_filename,
            ref_filename,
            remapper,
            remap_file=True,
        )

    def test_stereographic_array_to_latlon_array(self):
        """
        test horizontal interpolation from a stereographic grid to a
        destination lat/lon grid determined from config options 'lat' and
        'lon'.
        """

        weight_filename, out_filename, ref_filename = self.get_filenames(
            suffix='stereographic_to_latlon'
        )

        src_descriptor = self.get_stereographic_array_descriptor()
        dst_descriptor = self.get_latlon_array_descriptor()

        Lat = src_descriptor.coords['lat']['data']

        # now, let's make a more complicated field with more dimensions to
        # check
        in_field = np.reshape(Lat, (1, Lat.shape[0], Lat.shape[1], 1))
        in_field = in_field.repeat(3, axis=0)
        in_field = in_field.repeat(2, axis=3)

        dataset_dict = {
            'dims': ('dim0', 'y', 'x', 'dim3'),
            'coords': src_descriptor.coords,
            'data_vars': {
                'complicated': {
                    'dims': ('dim0', 'y', 'x', 'dim3'),
                    'data': in_field,
                }
            },
        }

        ds = xr.Dataset.from_dict(dataset_dict)
        in_filename = (
            f'{self.test_dir}/unmapped_stereographic_array_to_latlon_array.nc'
        )
        ds.to_netcdf(in_filename)

        remapper = self.build_remapper(
            src_descriptor, dst_descriptor, weight_filename
        )

        self.check_remap(
            in_filename, out_filename, ref_filename, remapper, remap_file=False
        )

    def test_latlon_file_to_point_collection(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        point collection.
        """

        weight_filename, out_filename, ref_filename = self.get_filenames(
            suffix='latlon_file_to_point_collection'
        )

        src_descriptor, latlon_grid_filename = (
            self.get_latlon_file_descriptor()
        )
        dst_descriptor = self.get_point_collection_descriptor()

        remapper = self.build_remapper(
            src_descriptor, dst_descriptor, weight_filename
        )
        self.check_remap(
            latlon_grid_filename,
            out_filename,
            ref_filename,
            remapper,
            remap_file=True,
        )
