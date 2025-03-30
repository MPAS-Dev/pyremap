# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
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

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        self.renormalizationThreshold = 0.01

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def get_mpas_cell_descriptor(self):
        mpasMeshFileName = str(self.datadir.join('mpasMesh.nc'))
        timeSeriesFileName = str(self.datadir.join('timeSeries.0002-01-01.nc'))
        descriptor = MpasCellMeshDescriptor(
            mpasMeshFileName, meshName='oQU240')

        return descriptor, mpasMeshFileName, timeSeriesFileName

    def get_mpas_edge_descriptor(self):
        mpasMeshFileName = str(self.datadir.join('mpasMesh.nc'))
        areaFileName = str(self.datadir.join('mpasAreaEdge.nc'))
        descriptor = MpasEdgeMeshDescriptor(
            mpasMeshFileName, meshName='oQU240')

        return descriptor, mpasMeshFileName, areaFileName

    def get_mpas_vertex_descriptor(self):
        mpasMeshFileName = str(self.datadir.join('mpasMesh.nc'))
        areaFileName = str(self.datadir.join('mpasAreaVertex.nc'))
        descriptor = MpasVertexMeshDescriptor(
            mpasMeshFileName, meshName='oQU240')

        return descriptor, mpasMeshFileName, areaFileName

    def get_latlon_file_descriptor(self):
        latLonGridFileName = str(self.datadir.join('SST_annual_1870-1900.nc'))
        descriptor = LatLonGridDescriptor.read(
            latLonGridFileName, latVarName='lat', lonVarName='lon')

        return descriptor, latLonGridFileName

    def get_latlon_array_descriptor(self):
        lat = np.linspace(-90., 90., 91)
        lon = np.linspace(-180., 180., 181)

        descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')
        return descriptor

    def get_latlon2d_file_descriptor(self):
        latLonGridFileName = str(self.datadir.join('SST_annual_1870-1900.nc'))
        ds = xr.open_dataset(latLonGridFileName)
        lon2d, lat2d = np.meshgrid(ds.lon.values, ds.lat.values)
        ds['lat2d'] = (('lat', 'lon'), lat2d)
        ds['lon2d'] = (('lat', 'lon'), lon2d)
        ds.lat2d.attrs['units'] = ds.lat.attrs['units']
        ds.lon2d.attrs['units'] = ds.lon.attrs['units']
        descriptor = LatLon2DGridDescriptor.read(
            ds=ds, latVarName='lat2d', lonVarName='lon2d')

        return descriptor, latLonGridFileName

    def get_point_collection_descriptor(self):
        mpasMeshFileName = str(self.datadir.join('mpasMesh.nc'))
        ds = xr.open_dataset(mpasMeshFileName)
        lats = ds.latCell.values
        lons = ds.lonCell.values

        descriptor = PointCollectionDescriptor(
            lats=lats,
            lons=lons,
            collectionName='mpasCellCenters',
            units='radians')
        return descriptor

    def get_stereographic_array_descriptor(self):

        # projection for BEDMAP2 and other common Antarctic data sets
        projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 '
                                 '+lon_0=0.0  +k_0=1.0 +x_0=0.0 +y_0=0.0 '
                                 '+ellps=WGS84')

        # a 61x51 cell map with 100 km resolution and
        xMax = 3000e3
        yMax = 2500e3
        res = 100e3
        nx = 2 * int(xMax / res) + 1
        ny = 2 * int(yMax / res) + 1
        x = np.linspace(-xMax, xMax, nx)
        y = np.linspace(-yMax, yMax, ny)
        meshName = f'{int(res * 1e-3)}km_Antarctic_stereo'
        descriptor = ProjectionGridDescriptor.create(
            projection, x, y, meshName)
        return descriptor

    def get_file_names(self, suffix):
        weightFileName = f'{self.test_dir}/weights_{suffix}.nc'
        outFileName = f'{self.test_dir}/remapped_{suffix}.nc'
        refFileName = f'{self.datadir}/ref_{suffix}.nc'
        return weightFileName, outFileName, refFileName

    def get_scrip_file_names(self, suffix):
        scripFileName = f'{self.test_dir}/scrip_{suffix}.nc'
        refFileName = f'{self.datadir}/ref_scrip_{suffix}.nc'
        return scripFileName, refFileName

    def build_remapper(self, sourceDescriptor, destinationDescriptor,
                       weightFileName):

        remapper = Remapper(sourceDescriptor, destinationDescriptor,
                            weightFileName)

        remapper.esmf_build_map(method='bilinear')

        assert os.path.exists(remapper.mappingFileName)

        return remapper

    def check_scrip(self, scripFileName, refFileName):

        dsRef = xr.open_dataset(refFileName)
        assert os.path.exists(scripFileName)
        dsScrip = xr.open_dataset(scripFileName)
        self.assertDatasetApproxEqual(dsScrip, dsRef)

    def check_remap(
            self,
            inFileName,
            outFileName,
            refFileName,
            remapper,
            remap_file=True):

        dropVars = [
            'lat_bnds',
            'lon_bnds',
            'gw',
            'area',
            'nvertices',
            'lat_vertices',
            'lon_vertices'
        ]
        dsRemappedFile = None
        if remap_file:
            # first, test interpolation with ncremap
            remapper.remap_file(
                inFileName=inFileName,
                outFileName=outFileName,
                replaceMpasFill=True)

            assert os.path.exists(outFileName)
            dsRemappedFile = xr.open_dataset(outFileName)
            # drop some extra vairables added by ncremap that aren't in the
            # reference data set
            dsRemappedFile = dsRemappedFile.drop_vars(
                [var for var in dropVars if var in dsRemappedFile]
            )
            newRefFilename = outFileName.replace('.nc', '_ref.nc')
            dsRemappedFile.to_netcdf(newRefFilename)

        # now, try in-memory remapping
        ds = xr.open_dataset(inFileName)
        dsRemapped = remapper.remap(ds, self.renormalizationThreshold)

        # see how we compare with the reference data set
        dsRef = xr.open_dataset(refFileName)
        dsRef = dsRef.drop_vars([var for var in dropVars if var in dsRef])
        if remap_file:
            self.assertDatasetApproxEqual(dsRemappedFile, dsRef)

        self.assertDatasetApproxEqual(dsRemapped, dsRef)

    def test_latlon_file_scrip(self):
        """
        test writing a SCRIP file for a lat/lon grid file
        """

        scripFileName, refFileName = \
            self.get_scrip_file_names(suffix='latlon_file')

        sourceDescriptor, _ = self.get_latlon_file_descriptor()
        print(f'Writing SCRIP file for lat/lon grid {scripFileName}')
        sourceDescriptor.to_scrip(scripFileName)

        self.check_scrip(scripFileName, refFileName)

    def test_latlon_array_scrip(self):
        """
        test writing a SCRIP file for a lat/lon grid array
        """

        scripFileName, refFileName = \
            self.get_scrip_file_names(suffix='latlon_array')

        sourceDescriptor = self.get_latlon_array_descriptor()
        print(f'Writing SCRIP file for lat/lon grid {scripFileName}')
        sourceDescriptor.to_scrip(scripFileName)

        self.check_scrip(scripFileName, refFileName)

    def test_latlon2d_scrip(self):
        """
        test writing a SCRIP file for a 2D lat/lon grid
        """

        scripFileName, refFileName = \
            self.get_scrip_file_names(suffix='latlon2d')

        sourceDescriptor, _ = self.get_latlon2d_file_descriptor()
        print(f'Writing SCRIP file for 2D lat/lon grid {scripFileName}')
        sourceDescriptor.to_scrip(scripFileName)

        self.check_scrip(scripFileName, refFileName)

    def test_mpas_cell_scrip(self):
        """
        test writing a SCRIP file for an MPAS cell mesh
        """

        scripFileName, refFileName = \
            self.get_scrip_file_names(suffix='mpas_cell')

        sourceDescriptor, _, _ = \
            self.get_mpas_cell_descriptor()
        print(f'Writing SCRIP file for MPAS cell mesh {scripFileName}')
        sourceDescriptor.to_scrip(scripFileName)

        self.check_scrip(scripFileName, refFileName)

    def test_mpas_edge_scrip(self):
        """
        test writing a SCRIP file for an MPAS edge mesh
        """

        scripFileName, refFileName = \
            self.get_scrip_file_names(suffix='mpas_edge')

        sourceDescriptor, _, _ = \
            self.get_mpas_edge_descriptor()
        print(f'Writing SCRIP file for MPAS edge mesh {scripFileName}')
        sourceDescriptor.to_scrip(scripFileName)

        self.check_scrip(scripFileName, refFileName)

    def test_mpas_vertex_scrip(self):
        """
        test writing a SCRIP file for an MPAS vertex mesh
        """

        scripFileName, refFileName = \
            self.get_scrip_file_names(suffix='mpas_vertex')

        sourceDescriptor, _, _ = \
            self.get_mpas_vertex_descriptor()
        print(f'Writing SCRIP file for MPAS vertex mesh {scripFileName}')
        sourceDescriptor.to_scrip(scripFileName)

        self.check_scrip(scripFileName, refFileName)

    def test_point_collection_scrip(self):
        """
        test writing a SCRIP file for a point collection
        """

        scripFileName, refFileName = \
            self.get_scrip_file_names(suffix='point_collection')

        sourceDescriptor = self.get_point_collection_descriptor()
        print(f'Writing SCRIP file for point collection {scripFileName}')
        sourceDescriptor.to_scrip(scripFileName)

        self.check_scrip(scripFileName, refFileName)

    def test_stereographic_scrip(self):
        """
        test writing a SCRIP file for a stereographic projection grid
        """

        scripFileName, refFileName = \
            self.get_scrip_file_names(suffix='stereographic')

        sourceDescriptor = self.get_stereographic_array_descriptor()
        print(f'Writing SCRIP file for stereographic grid {scripFileName}')
        sourceDescriptor.to_scrip(scripFileName)

        self.check_scrip(scripFileName, refFileName)

    def test_mpas_cell_to_latlon(self):
        """
        test horizontal interpolation from an MPAS cell mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords
        """

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='mpas_cell_to_latlon')

        sourceDescriptor, _, timeSeriesFileName = \
            self.get_mpas_cell_descriptor()
        destinationDescriptor, _ = self.get_latlon_file_descriptor()

        remapper = self.build_remapper(
            sourceDescriptor, destinationDescriptor, weightFileName)
        self.check_remap(
            timeSeriesFileName,
            outFileName,
            refFileName,
            remapper,
            remap_file=True
        )

    def test_mpas_edge_to_latlon(self):
        """
        test horizontal interpolation from an MPAS edge mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords
        """

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='mpas_edge_to_latlon')

        sourceDescriptor, _, areaFileName = \
            self.get_mpas_edge_descriptor()
        destinationDescriptor, _ = self.get_latlon_file_descriptor()

        remapper = self.build_remapper(
            sourceDescriptor, destinationDescriptor, weightFileName)
        self.check_remap(
            areaFileName,
            outFileName,
            refFileName,
            remapper,
            remap_file=True
        )

    def test_mpas_vertex_to_latlon(self):
        """
        test horizontal interpolation from an MPAS vertex mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords
        """

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='mpas_vertex_to_latlon')

        sourceDescriptor, _, areaFileName = \
            self.get_mpas_vertex_descriptor()
        destinationDescriptor, _ = self.get_latlon_file_descriptor()

        remapper = self.build_remapper(
            sourceDescriptor, destinationDescriptor, weightFileName)
        self.check_remap(
            areaFileName,
            outFileName,
            refFileName,
            remapper,
            remap_file=True
        )

    def test_latlon_file_to_latlon_array(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        lat/lon grid determined from config options 'lat' and 'lon'.
        """

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='latlon_file_to_latlon_array')

        sourceDescriptor, latLonGridFileName = \
            self.get_latlon_file_descriptor()
        destinationDescriptor = self.get_latlon_array_descriptor()

        remapper = self.build_remapper(
            sourceDescriptor, destinationDescriptor, weightFileName)
        self.check_remap(
            latLonGridFileName,
            outFileName,
            refFileName,
            remapper,
            remap_file=True
        )

    def test_mpas_cell_to_stereographic(self):
        """
        test horizontal interpolation from an MPAS cell mesh to a destination
        stereographic grid.
        """

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='mpas_cell_to_stereographic')

        sourceDescriptor, _, timeSeriesFileName = \
            self.get_mpas_cell_descriptor()
        destinationDescriptor = self.get_stereographic_array_descriptor()

        remapper = self.build_remapper(
            sourceDescriptor, destinationDescriptor, weightFileName)

        # ncremap doesn't support stereographic grids so just check the
        # Remapper
        self.check_remap(
            timeSeriesFileName,
            outFileName,
            refFileName,
            remapper,
            remap_file=True
        )

    def test_latlon_to_stereographic(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        stereographic grid.
        """

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='latlon_to_stereographic')

        sourceDescriptor, latLonGridFileName = \
            self.get_latlon_file_descriptor()
        destinationDescriptor = self.get_stereographic_array_descriptor()

        remapper = self.build_remapper(
            sourceDescriptor, destinationDescriptor, weightFileName)

        # ncremap doesn't support stereographic grids so just check the
        # Remapper
        self.check_remap(
            latLonGridFileName,
            outFileName,
            refFileName,
            remapper,
            remap_file=True)

    def test_stereographic_array_to_latlon_array(self):
        """
        test horizontal interpolation from a stereographic grid to a
        destination lat/lon grid determined from config options 'lat' and
        'lon'.
        """

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='stereographic_to_latlon')

        sourceDescriptor = self.get_stereographic_array_descriptor()
        destinationDescriptor = self.get_latlon_array_descriptor()

        Lat = sourceDescriptor.coords['lat']['data']

        # now, let's make a more complicated field with more dimensions to
        # check
        inField = np.reshape(Lat, (1, Lat.shape[0], Lat.shape[1], 1))
        inField = inField.repeat(3, axis=0)
        inField = inField.repeat(2, axis=3)

        datasetDict = {'dims': ('dim0', 'y', 'x', 'dim3'),
                       'coords': sourceDescriptor.coords,
                       'data_vars': {'complicated':
                                     {'dims': ('dim0', 'y', 'x', 'dim3'),
                                      'data': inField}}}

        ds = xr.Dataset.from_dict(datasetDict)
        inFileName = (
            f'{self.test_dir}/unmapped_stereographic_array_to_latlon_'
            f'array.nc'
        )
        ds.to_netcdf(inFileName)

        remapper = self.build_remapper(sourceDescriptor, destinationDescriptor,
                                       weightFileName)

        self.check_remap(inFileName, outFileName, refFileName,
                         remapper, remap_file=False)

    def test_latlon_file_to_point_collection(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        point collection.
        """

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='latlon_file_to_point_collection')

        sourceDescriptor, latLonGridFileName = \
            self.get_latlon_file_descriptor()
        destinationDescriptor = self.get_point_collection_descriptor()

        remapper = self.build_remapper(
            sourceDescriptor, destinationDescriptor, weightFileName)
        self.check_remap(
            latLonGridFileName,
            outFileName,
            refFileName,
            remapper,
            remap_file=True
        )
