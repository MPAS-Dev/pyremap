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

import json
import os
import subprocess
import sys
import warnings
from distutils.spawn import find_executable
from subprocess import check_output
from tempfile import TemporaryDirectory

import numpy
import xarray as xr
from scipy.sparse import csr_matrix

from pyremap.descriptor import (
    LatLon2DGridDescriptor,
    LatLonGridDescriptor,
    MpasEdgeMeshDescriptor,
    MpasMeshDescriptor,
    PointCollectionDescriptor,
    ProjectionGridDescriptor,
)


class Remapper(object):
    """
    A class for remapping fields using a given mapping file.  The weights and
    indices from the mapping file can be loaded once and reused multiple times
    to map several fields between the same source and destination grids.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, sourceDescriptor, destinationDescriptor,
                 mappingFileName=None):
        """
        Create the remapper and read weights and indices from the given file
        for later used in remapping fields.

        Parameters
        ----------
        sourceDescriptor : ``shared.grid.MeshDescriptor``
            An object used to write a scrip file and to determine the type of
            the source mesh or grid.

        destinationDescriptor : ``shared.grid.MeshDescriptor``
            An object used to write a scrip files and to determine the type of
            the destination mesh or grid.

        mappingFileName : str, optional
            The path where the mapping file containing interpolation weights
            and indices will be written and/or read.  If ``None``,
            no interpolation is performed and data sets are returned unchanged.
            This is useful if the source and destination grids are determined
            to be the same (though the Remapper does not attempt to determine
            if this is the case).
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if isinstance(sourceDescriptor, PointCollectionDescriptor):
            raise TypeError("sourceDescriptor of type "
                            "PointCollectionDescriptor is not supported.")
        if not isinstance(sourceDescriptor,
                          (MpasMeshDescriptor, MpasEdgeMeshDescriptor,
                           LatLonGridDescriptor,
                           LatLon2DGridDescriptor, ProjectionGridDescriptor)):
            raise TypeError("sourceDescriptor is not of a recognized type.")

        if not isinstance(destinationDescriptor,
                          (MpasMeshDescriptor, MpasEdgeMeshDescriptor,
                           LatLonGridDescriptor,
                           LatLon2DGridDescriptor, ProjectionGridDescriptor,
                           PointCollectionDescriptor)):
            raise TypeError(
                "destinationDescriptor is not of a recognized type.")

        self.sourceDescriptor = sourceDescriptor
        self.destinationDescriptor = destinationDescriptor
        self.mappingFileName = mappingFileName

        self.mappingLoaded = False

    def build_mapping_file(self, method='bilinear',  # noqa: C901
                           additionalArgs=None, logger=None, mpiTasks=1,
                           tempdir=None, esmf_path=None,
                           esmf_parallel_exec=None, extrap_method=None):
        """
        Given a source file defining either an MPAS mesh or a lat-lon grid and
        a destination file or set of arrays defining a lat-lon grid, constructs
        a mapping file used for interpolation between the source and
        destination grids.

        Parameters
        ----------
        method : {'bilinear', 'neareststod', 'conserve'}, optional
            The method of interpolation used, see documentation for
            `ESMF_RegridWeightGen` for details.

        additionalArgs : list of str, optional
            A list of additional arguments to ``ESMF_RegridWeightGen``

        logger : ``logging.Logger``, optional
            A logger to which ncclimo output should be redirected

        mpiTasks : int, optional
            The number of MPI tasks (a number > 1 implies that
            ESMF_RegridWeightGen will be called with ``mpirun``)

        tempdir : str, optional
            A temporary directory.  By default, a temporary directory is
            created, typically in ``/tmp`` but on some systems such as compute
            nodes this may not be visible to all processors in the subsequent
            ``ESMF_RegridWeightGen`` call

        esmf_path : str, optional
            A path to a system build of ESMF (containing a 'bin' directory with
            the ESMF tools).  By default, ESMF tools are found in the conda
            environment

        esmf_parallel_exec : {'srun', 'mpirun}, optional
            The name of the parallel executable to use to launch ESMF tools.
            By default, 'mpirun' from the conda environment is used

        extrap_method : {'neareststod', 'nearestidavg','creep'}, optional
            The method used to extrapolate unmapped destination locations

        Raises
        ------
        OSError
            If ``ESMF_RegridWeightGen`` is not in the system path.

        ValueError
            If sourceDescriptor or destinationDescriptor is of an unknown type
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if isinstance(self.destinationDescriptor,
                      PointCollectionDescriptor) and \
                method not in ['bilinear', 'neareststod']:
            raise ValueError("method {} not supported for destination "
                             "grid of type PointCollectionDescriptor."
                             "".format(method))

        if self.mappingFileName is None or \
                os.path.exists(self.mappingFileName):
            # a valid weight file already exists, so nothing to do
            return

        if esmf_path is not None:
            # use the system build of ESMF
            rwgPath = os.path.join(esmf_path, 'bin', 'ESMF_RegridWeightGen')

        else:
            rwgPath = find_executable('ESMF_RegridWeightGen')

            if rwgPath is None:
                raise OSError('ESMF_RegridWeightGen not found. Make sure esmf '
                              'package is installed: \n'
                              'conda install esmf\n'
                              'Note: this presumes use of the conda-forge '
                              'channel.')

        # Write source and destination SCRIP files in temporary locations
        if tempdir is None:
            tempobj = TemporaryDirectory()
            tempdir = tempobj.name
        else:
            tempobj = None

        sourceFileName = '{}/src_mesh.nc'.format(tempdir)
        destinationFileName = '{}/dst_mesh.nc'.format(tempdir)

        src_loc = 'center'
        src_file_format = 'scrip'
        if isinstance(self.sourceDescriptor,
                      (MpasMeshDescriptor, MpasEdgeMeshDescriptor)):
            src_file_format = 'esmf'

        if isinstance(self.sourceDescriptor, MpasMeshDescriptor) and \
                self.sourceDescriptor.vertices:
            if 'conserve' in method:
                raise ValueError('Can\'t remap from MPAS vertices with '
                                 'conservative methods')
            src_loc = 'corner'

        dst_loc = 'center'
        dst_file_format = 'scrip'
        if isinstance(self.destinationDescriptor,
                      (MpasMeshDescriptor, MpasEdgeMeshDescriptor)):
            dst_file_format = 'esmf'

        if isinstance(self.destinationDescriptor, MpasMeshDescriptor) and \
                self.destinationDescriptor.vertices:
            if 'conserve' in method:
                raise ValueError('Can\'t remap to MPAS vertices with '
                                 'conservative methods')
            dst_loc = 'corner'

        if src_file_format == 'scrip':
            self.sourceDescriptor.to_scrip(sourceFileName)
        elif src_file_format == 'esmf':
            self.sourceDescriptor.to_esmf(sourceFileName)
        else:
            raise ValueError('Unexpected file format {}'.format(
                src_file_format))

        if dst_file_format == 'scrip':
            self.destinationDescriptor.to_scrip(destinationFileName)
        elif dst_file_format == 'esmf':
            self.destinationDescriptor.to_esmf(destinationFileName)
        else:
            raise ValueError('Unexpected file format {}'.format(
                dst_file_format))

        args = [rwgPath,
                '--source', sourceFileName,
                '--destination', destinationFileName,
                '--weight', self.mappingFileName,
                '--method', method,
                '--netcdf4',
                '--no_log']

        if extrap_method is not None:
            args.extend(['--extrap_method', extrap_method])

        if src_file_format == 'esmf':
            args.extend(['--src_loc', src_loc])
        if dst_file_format == 'esmf':
            args.extend(['--dst_loc', dst_loc])

        parallel_args = []

        if esmf_parallel_exec is not None:
            # use the specified parallel executable
            parallel_args = esmf_parallel_exec.split(' ')

            if 'srun' in esmf_parallel_exec:
                parallel_args.extend(['-n', '{}'.format(mpiTasks)])
            else:
                # presume mpirun syntax
                parallel_args.extend(['-np', '{}'.format(mpiTasks)])

        elif 'CONDA_PREFIX' in os.environ and mpiTasks > 1:
            # this is a conda environment, so we need to find out if esmf
            # needs mpirun or not
            conda_args = ['conda', 'list', 'esmf', '--json']
            output = check_output(conda_args).decode("utf-8")
            output = json.loads(output)
            build_string = output[0]['build_string']

            if 'mpi_mpich' in build_string or 'mpi_openmpi' in build_string:
                # esmf was installed with MPI, so we should use mpirun
                mpirun_path = '{}/bin/mpirun'.format(
                    os.environ['CONDA_PREFIX'])
                parallel_args = [mpirun_path, '-np', '{}'.format(mpiTasks)]
            else:
                # esmf was installed without MPI, so we shouldn't try to
                # use it
                warnings.warn('Requesting {} MPI tasks but the MPI version'
                              ' of ESMF is not installed'.format(mpiTasks))

        args = parallel_args + args

        if self.sourceDescriptor.regional:
            args.append('--src_regional')

        if self.destinationDescriptor.regional:
            args.append('--dst_regional')

        if self.sourceDescriptor.regional or \
                self.destinationDescriptor.regional:
            args.append('--ignore_unmapped')

        if additionalArgs is not None:
            args.extend(additionalArgs)

        if logger is None:
            _print_running(args, fn=print)
            # make sure any output is flushed before we add output from the
            # subprocess
            sys.stdout.flush()
            sys.stderr.flush()

            # throw out the standard output from ESMF_RegridWeightGen, as it's
            # rather verbose but keep stderr
            with open(os.devnull, 'wb') as DEVNULL:
                subprocess.check_call(args, stdout=DEVNULL)

        else:
            _print_running(args, fn=logger.info)
            for handler in logger.handlers:
                handler.flush()

            process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            # throw out the standard output from ESMF_RegridWeightGen, as it's
            # rather verbose but keep stderr
            if stderr:
                stderr = stderr.decode('utf-8')
                for line in stderr.split('\n'):
                    logger.error(line)

            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode,
                                                    ' '.join(args))

        if tempobj is not None:
            tempobj.cleanup()

    def remap_file(self, inFileName, outFileName,  # noqa: C901
                   variableList=None, overwrite=False, renormalize=None,
                   logger=None, replaceMpasFill=False, parallel_exec=None):
        """
        Given a source file defining either an MPAS mesh or a lat-lon grid and
        a destination file or set of arrays defining a lat-lon grid, constructs
        a mapping file used for interpolation between the source and
        destination grids.

        Parameters
        ----------
        inFileName : str
            The path to the file containing a data set on the source grid

        outFileName : str
            The path where the data on the destination grid should be written

        variableList : list of str, optional
            A list of variables to be mapped.  By default, all variables are
            mapped

        overwrite : bool, optional
            Whether the destination file should be overwritten if it already
            exists. If `False`, and the destination file is already present,
            the function does nothing and returns immediately

        renormalize : float, optional
            A threshold to use to renormalize the data

        logger : ``logging.Logger``, optional
            A logger to which ncclimo output should be redirected

        replaceMpasFill : bool, optional
            For MPAS meshes, whether add a ``_FillValue`` attribute (missing
            from MPAS output).  If this has been handled before the call,
            replacing the fill value again may cause errors.

        parallel_exec : {'srun'}, optional
            The name of the parallel executable to use to launch ncremap.
            By default, none is used.

        Raises
        ------
        OSError
            If ``ncremap`` is not in the system path.

        ValueError
            If ``mappingFileName`` is ``None`` (meaning no remapping is
            needed).
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if self.mappingFileName is None:
            raise ValueError('No mapping file was given because remapping is '
                             'not necessary. The calling\n'
                             'code should simply use the constents of {} '
                             'directly.'.format(inFileName))

        if not overwrite and os.path.exists(outFileName):
            # a remapped file already exists, so nothing to do
            return

        if isinstance(self.sourceDescriptor, PointCollectionDescriptor):
            raise TypeError('Source grid is a point collection, which is not'
                            'supported.')

        if find_executable('ncremap') is None:
            raise OSError('ncremap not found. Make sure the latest nco '
                          'package is installed: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        if parallel_exec is not None:
            # use the specified parallel executable
            args = parallel_exec.split(' ')
        else:
            args = list()

        args.extend(['ncremap',
                     '-m', self.mappingFileName,
                     '--vrb=1'])

        regridArgs = []

        if isinstance(self.sourceDescriptor, MpasMeshDescriptor):
            if self.sourceDescriptor.vertices:
                regridArgs.extend(['--rgr col_nm=nVertices'])
            else:
                args.extend(['-P', 'mpas'])
                if not replaceMpasFill:
                    # the -C (climatology) flag prevents ncremap from trying to
                    # add a _FillValue attribute that might already be present
                    # and quits with an error
                    args.append('-C')

        if isinstance(self.sourceDescriptor, MpasEdgeMeshDescriptor):
            regridArgs.extend(['--rgr col_nm=nEdges'])

        if isinstance(self.sourceDescriptor,
                      (MpasMeshDescriptor, MpasEdgeMeshDescriptor)) and \
                renormalize is not None:
            # we also want to make sure cells that receive no data are
            # marked with fill values, even if the source MPAS data
            # doesn't have a fill value
            args.append('--add_fill_value')

        if variableList is not None:
            args.extend(['-v', ','.join(variableList)])

        if renormalize is not None:
            regridArgs.append('--renormalize={}'.format(renormalize))

        if isinstance(self.sourceDescriptor, LatLonGridDescriptor):
            regridArgs.extend(['--rgr lat_nm={}'.format(
                self.sourceDescriptor.latVarName),
                '--rgr lon_nm={}'.format(
                self.sourceDescriptor.lonVarName)])
        elif isinstance(self.sourceDescriptor, ProjectionGridDescriptor):
            regridArgs.extend(['--rgr lat_nm={}'.format(
                self.sourceDescriptor.yVarName),
                '--rgr lon_nm={}'.format(
                self.sourceDescriptor.xVarName)])

        if isinstance(self.destinationDescriptor, LatLonGridDescriptor):
            regridArgs.extend([
                f'--rgr lat_nm_out={self.destinationDescriptor.latVarName}',
                f'--rgr lon_nm_out={self.destinationDescriptor.lonVarName}',
                f'--rgr lat_dmn_nm={self.destinationDescriptor.latVarName}',
                f'--rgr lon_dmn_nm={self.destinationDescriptor.lonVarName}'])
        elif isinstance(self.destinationDescriptor, ProjectionGridDescriptor):
            regridArgs.extend(['--rgr lat_dmn_nm={}'.format(
                self.destinationDescriptor.yVarName),
                '--rgr lon_dmn_nm={}'.format(
                self.destinationDescriptor.xVarName),
                '--rgr lat_nm_out=lat', '--rgr lon_nm_out=lon'])
        if isinstance(self.destinationDescriptor, PointCollectionDescriptor):
            regridArgs.extend(['--rgr lat_nm_out=lat', '--rgr lon_nm_out=lon'])

        if len(regridArgs) > 0:
            args.extend(['-R', ' '.join(regridArgs)])

        # set an environment variable to make sure we're not using czender's
        # local version of NCO instead of one we have intentionally loaded
        env = os.environ.copy()
        env['NCO_PATH_OVERRIDE'] = 'No'

        args.extend([inFileName, outFileName])

        if logger is None:
            # make sure any output is flushed before we add output from the
            # subprocess
            sys.stdout.flush()
            sys.stderr.flush()

            _print_running(args, fn=print)
            subprocess.check_call(args, env=env)
        else:
            _print_running(args, fn=logger.info)
            for handler in logger.handlers:
                handler.flush()

            process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, env=env)
            stdout, stderr = process.communicate()

            if stdout:
                stdout = stdout.decode('utf-8')
                for line in stdout.split('\n'):
                    logger.info(line)
            if stderr:
                stderr = stderr.decode('utf-8')
                for line in stderr.split('\n'):
                    logger.error(line)

            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode,
                                                    ' '.join(args))

    def remap(self, ds, renormalizationThreshold=None):
        """
        Given a source data set, returns a remapped version of the data set,
        possibly masked and renormalized.

        Parameters
        ----------
        ds : ``xarray.Dataset`` or ``xarray.DataArray``
            The dimention(s) along ``self.sourceDimNames`` must match
            ``self.src_grid_dims`` read from the mapping file.

        renormalizationThreshold : float, optional
            The minimum weight of a denstination cell after remapping, below
            which it is masked out, or ``None`` for no renormalization and
            masking.

        Returns
        -------
        remappedDs : `xarray.Dataset`` or ``xarray.DataArray``
            Returns a remapped data set (or data array) where dimensions other
            than ``self.sourceDimNames`` are the same as in ``ds`` and the
            dimension(s) given by ``self.sourceDimNames`` have been replaced by
            ``self.destinationDimNames``.

        Raises
        ------
        ValueError
            If the size of ``self.sourceDimNames`` in ``ds`` do not match the
            source dimensions read in from the mapping file
            (``self.src_grid_dims``).
        TypeError
            If ds is not an ``xarray.Dataset`` or ``xarray.DataArray`` object
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if self.mappingFileName is None:
            # No remapping is needed
            return ds

        self._load_mapping()

        for index, dim in enumerate(self.sourceDescriptor.dims):
            if self.src_grid_dims[index] != ds.sizes[dim]:
                raise ValueError('data set and remapping source dimension {} '
                                 'don\'t have the same size: {} != {}'.format(
                                     dim, self.src_grid_dims[index],
                                     ds.sizes[dim]))

        if isinstance(ds, xr.DataArray):
            remappedDs = self._remap_data_array(ds, renormalizationThreshold)
        elif isinstance(ds, xr.Dataset):
            drop = []
            for var in ds.data_vars:
                if self._check_drop(ds[var]):
                    drop.append(var)
            remappedDs = ds.drop_vars(drop)
            remappedDs = remappedDs.map(self._remap_data_array,
                                        keep_attrs=True,
                                        args=(renormalizationThreshold,))
        else:
            raise TypeError('ds not an xarray Dataset or DataArray.')

        # Update history attribute of netCDF file
        if 'history' in remappedDs.attrs:
            newhist = '\n'.join([remappedDs.attrs['history'],
                                 ' '.join(sys.argv[:])])
        else:
            newhist = sys.argv[:]
        remappedDs.attrs['history'] = newhist

        remappedDs.attrs['meshName'] = self.destinationDescriptor.meshName

        return remappedDs

    def _load_mapping(self):
        """
        Load weights and indices from a mapping file, if this has not already
        been done
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if self.mappingLoaded:
            return

        dsMapping = xr.open_dataset(self.mappingFileName)
        n_a = dsMapping.dims['n_a']
        n_b = dsMapping.dims['n_b']

        nSourceDims = len(self.sourceDescriptor.dims)
        src_grid_rank = dsMapping.dims['src_grid_rank']
        nDestinationDims = len(self.destinationDescriptor.dims)
        dst_grid_rank = dsMapping.dims['dst_grid_rank']

        # check that the mapping file has the right number of dimensions
        if nSourceDims != src_grid_rank or \
                nDestinationDims != dst_grid_rank:
            raise ValueError('The number of source and/or '
                             'destination dimensions does not\n'
                             'match the expected number of source and '
                             'destination dimensions in the mapping\n'
                             'file. {} != {} and/or {} != {}'.format(
                                 nSourceDims, src_grid_rank,
                                 nDestinationDims, dst_grid_rank))

        # grid dimensions need to be reversed because they are in Fortran order
        self.src_grid_dims = dsMapping['src_grid_dims'].values[::-1]
        self.dst_grid_dims = dsMapping['dst_grid_dims'].values[::-1]

        # now, check that each source and destination dimension is right
        for index in range(len(self.sourceDescriptor.dims)):
            dim = self.sourceDescriptor.dims[index]
            dimSize = self.sourceDescriptor.dimSize[index]
            checkDimSize = self.src_grid_dims[index]
            if dimSize != checkDimSize:
                raise ValueError('source mesh descriptor and remapping source '
                                 'dimension {} don\'t have the same size: \n'
                                 '{} != {}'.format(dim, dimSize, checkDimSize))
        for index in range(len(self.destinationDescriptor.dims)):
            dim = self.destinationDescriptor.dims[index]
            dimSize = self.destinationDescriptor.dimSize[index]
            checkDimSize = self.dst_grid_dims[index]
            if dimSize != checkDimSize:
                raise ValueError('dest. mesh descriptor and remapping dest. '
                                 'dimension {} don\'t have the same size: \n'
                                 '{} != {}'.format(dim, dimSize, checkDimSize))

        self.frac_b = dsMapping['frac_b'].values

        col = dsMapping['col'].values - 1
        row = dsMapping['row'].values - 1
        S = dsMapping['S'].values
        self.matrix = csr_matrix((S, (row, col)), shape=(n_b, n_a))

        self.mappingLoaded = True

    def _check_drop(self, dataArray):
        sourceDims = self.sourceDescriptor.dims

        sourceDimsInArray = [dim in dataArray.dims for dim in sourceDims]

        return (numpy.any(sourceDimsInArray) and not
                numpy.all(sourceDimsInArray))

    def _remap_data_array(self, dataArray, renormalizationThreshold):
        """
        Remap a single xarray data array
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        sourceDims = self.sourceDescriptor.dims
        destDims = self.destinationDescriptor.dims

        sourceDimsInArray = [dim in dataArray.dims for dim in sourceDims]

        if not numpy.any(sourceDimsInArray):
            # no remapping is needed
            return dataArray

        if not numpy.all(sourceDimsInArray):
            # no remapping is possible so the variable array should have been
            # dropped
            raise ValueError('Data array with some (but not all) required '
                             'source dims cannot be remapped\n'
                             'and should have been dropped.')

        # make a list of dims and remapAxes
        dims = []
        remapAxes = []
        destDimsAdded = False
        for index, dim in enumerate(dataArray.dims):
            if dim in sourceDims:
                remapAxes.append(index)
                if not destDimsAdded:
                    dims.extend(destDims)
                    destDimsAdded = True
            else:
                dims.append(dim)

        # make a dict of coords
        coordDict = {}
        # copy unmodified coords
        for coord in dataArray.coords:
            sourceDimInCoord = numpy.any([dim in dataArray.coords[coord].dims
                                          for dim in sourceDims])
            if not sourceDimInCoord:
                coordDict[coord] = {'dims': dataArray.coords[coord].dims,
                                    'data': dataArray.coords[coord].values}

        # add dest coords
        coordDict.update(self.destinationDescriptor.coords)

        # remap the values
        field = dataArray.values
        mask = numpy.isnan(field)
        if numpy.count_nonzero(mask) > 0:
            field = numpy.ma.masked_array(field, mask)
        remappedField = self._remap_numpy_array(field, remapAxes,
                                                renormalizationThreshold)

        arrayDict = {'coords': coordDict,
                     'attrs': dataArray.attrs,
                     'dims': dims,
                     'data': remappedField,
                     'name': dataArray.name}

        # make a new data array
        remappedArray = xr.DataArray.from_dict(arrayDict)

        return remappedArray

    def _remap_numpy_array(self, inField, remapAxes,
                           renormalizationThreshold):
        """
        Remap a single numpy array
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # permute the dimensions of inField so the axes to remap are first,
        # then flatten the remapping and the extra dimensions separately for
        # the matrix multiply
        extraAxes = [axis for axis in numpy.arange(inField.ndim)
                     if axis not in remapAxes]

        newShape = [numpy.prod([inField.shape[axis] for axis in remapAxes])]
        if len(extraAxes) > 0:
            extraShape = [inField.shape[axis] for axis in extraAxes]
            newShape.append(numpy.prod(extraShape))
        else:
            extraShape = []
            newShape.append(1)

        permutedAxes = remapAxes + extraAxes

        # permute axes so the remapped dimension(s) come first and "flatten"
        # the remapping dimension
        inField = inField.transpose(permutedAxes).reshape(newShape)

        masked = (isinstance(inField, numpy.ma.MaskedArray) and
                  renormalizationThreshold is not None)
        if masked:
            inMask = numpy.array(numpy.logical_not(inField.mask), float)
            outField = self.matrix.dot(inMask * inField)
            outMask = self.matrix.dot(inMask)
            mask = outMask > renormalizationThreshold
        else:
            outField = self.matrix.dot(inField)
            # make frac_b match the shape of outField
            outMask = numpy.reshape(self.frac_b, (len(self.frac_b), 1)).repeat(
                newShape[1], axis=1)
            mask = outMask > 0.

        # normalize the result based on outMask
        outField[mask] /= outMask[mask]
        outField = numpy.ma.masked_array(outField,
                                         mask=numpy.logical_not(mask))

        destRemapDimCount = len(self.dst_grid_dims)
        outDimCount = len(extraShape) + destRemapDimCount

        # "unflatten" the remapped dimension(s)
        destShape = list(self.dst_grid_dims) + extraShape
        outField = numpy.reshape(outField, destShape)

        # "unpermute" the axes to be in the expected order
        index = numpy.amin(remapAxes)
        unpermuteAxes = list(numpy.arange(destRemapDimCount, outDimCount))
        unpermuteAxes = (unpermuteAxes[0:index] +
                         list(numpy.arange(destRemapDimCount)) +
                         unpermuteAxes[index:])
        outField = numpy.transpose(outField, axes=unpermuteAxes)

        return outField


def _print_running(args, fn):
    print_args = []
    for arg in args:
        if ' ' in arg:
            arg = '"{}"'.format(arg)
        print_args.append(arg)
    fn('running: {}'.format(' '.join(print_args)))
