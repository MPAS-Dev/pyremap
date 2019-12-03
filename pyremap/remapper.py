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

import subprocess
import tempfile
import os
from distutils.spawn import find_executable
import numpy
from scipy.sparse import csr_matrix
import xarray
import sys

from pyremap import PointCollectionDescriptor, MeshDescriptor, write_netcdf


class Remapper(object):
    """
    A class for remapping fields using a given mapping file.  The weights and
    indices from the mapping file can be loaded once and reused multiple times
    to map several fields between the same source and destination grids.

    Attributes
    ----------
    src_descrip : ``shared.grid.MeshDescriptor``
        An object used to write a scrip file and to determine the type of
        the source mesh or grid.

    dst_descrip : ``shared.grid.MeshDescriptor``
        An object used to write a scrip files and to determine the type of
        the destination mesh or grid.

    mapping_filename : str
        The path where the mapping file containing interpolation weights
        and indices will be written and/or read.  If ``None``,
        no interpolation is performed and data sets are returned unchanged.
        This is useful if the source and destination grids are determined
        to be the same (though the Remapper does not attempt to determine
        if this is the case).
    """

    def __init__(self, src_descrip, dst_descrip,
                 mapping_filename=None):
        """
        Create the remapper and read weights and indices from the given file
        for later used in remapping fields.

        Parameters
        ----------
        src_descrip : MeshDescriptor
            An object used to write a scrip file and to determine the type of
            the source mesh or grid.

        dst_descrip : MeshDescriptor
            An object used to write a scrip files and to determine the type of
            the destination mesh or grid.

        mapping_filename : str, optional
            The path where the mapping file containing interpolation weights
            and indices will be written and/or read.  If ``None``,
            no interpolation is performed and data sets are returned unchanged.
            This is useful if the source and destination grids are determined
            to be the same (though ``Remapper`` does not attempt to determine
            if this is the case).
        """

        if isinstance(src_descrip, PointCollectionDescriptor):
            raise TypeError("PointCollectionDescriptor can only be a "
                            "destination.")

        self.src_descrip = src_descrip
        self.dst_descrip = dst_descrip
        self.mapping_filename = mapping_filename

        self._mapping_loaded = False

    def build_mapping_file(self, method='bilinear', software='esmf',
                           additional_args=None, logger=None, mpitasks=1):
        """
        Given a source file defining either an MPAS mesh or a lat-lon grid and
        a destination file or set of arrays defining a lat-lon grid, constructs
        a mapping file used for interpolation between the source and
        destination grids.

        Parameters
        ----------
        method : {'bilinear', 'neareststod', 'conserve'}, optional
            The method of interpolation used, see documentation for
            ``ESMF_RegridWeightGen`` for details.

        software : {'esmf'}, optional
            The software for creating the mapping file, currently only
            ``ESMF_RegridWeightGen`` is supported.

        additional_args : list of str, optional
            A list of additional arguments to ``ESMF_RegridWeightGen``

        logger : logging.Logger, optional
            A logger to which ncclimo output should be redirected

        mpitasks : int, optional
            The number of MPI tasks (a number > 1 implies that
            ESMF_RegridWeightGen will be called with ``mpirun``)
        """

        if isinstance(self.dst_descrip,
                      PointCollectionDescriptor) and \
                method not in ['bilinear', 'neareststod']:
            raise ValueError("method {} not supported for destination "
                             "grid of type PointCollectionDescriptor."
                             "".format(method))

        if self.mapping_filename is None or \
                os.path.exists(self.mapping_filename):
            # a valid weight file already exists, so nothing to do
            return

        if software != 'esmf':
            raise ValueError('Unknown software {}'.format(software))
        erwg_path = find_executable('ESMF_RegridWeightGen')

        if erwg_path is None:
            raise OSError('ESMF_RegridWeightGen not found. Make sure esmf '
                          'package is installed via\n'
                          'latest nco: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        # Write source and destination SCRIP files in temporary locations
        src_filename = _get_temp_path()
        dst_filename = _get_temp_path()
        self.src_descrip.to_scrip(src_filename)
        self.dst_descrip.to_scrip(dst_filename)

        args = [erwg_path,
                '--source', src_filename,
                '--destination', dst_filename,
                '--weight', self.mapping_filename,
                '--method', method,
                '--netcdf4',
                '--no_log']

        if mpitasks > 1:
            args = ['mpirun', '-np', '{}'.format(mpitasks)] + args

        if self.src_descrip.regional:
            args.append('--src_regional')

        if self.dst_descrip.regional:
            args.append('--dst_regional')

        if self.src_descrip.regional or \
                self.dst_descrip.regional:
            args.append('--ignore_unmapped')

        if additional_args is not None:
            args.extend(additional_args)

        if logger is None:
            print('running: {}'.format(' '.join(args)))
            # make sure any output is flushed before we add output from the
            # subprocess
            sys.stdout.flush()
            sys.stderr.flush()

            # throw out the standard output from ESMF_RegridWeightGen, as it's
            # rather verbose but keep stderr
            DEVNULL = open(os.devnull, 'wb')
            subprocess.check_call(args, stdout=DEVNULL)

        else:
            logger.info('running: {}'.format(' '.join(args)))
            for handler in logger.handlers:
                handler.flush()

            process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            # throw out the standard output from ESMF_RegridWeightGen, as it's
            # rather verbose but keep stderr
            if stderr:
                for line in stderr.split('\n'):
                    logger.error(line)

            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode,
                                                    ' '.join(args))

        # remove the temporary SCRIP files
        os.remove(src_filename)
        os.remove(dst_filename)

    def remap_file(self, infilename, outfilename, varlist=None,
                   overwrite=False, renormalize=None, logger=None):
        """
        Given a source file defining either an MPAS mesh or a lat-lon grid and
        a destination file or set of arrays defining a lat-lon grid, constructs
        a mapping file used for interpolation between the source and
        destination grids.

        Parameters
        ----------
        infilename : str
            The path to the file containing a data set on the source grid

        outfilename : str
            The path where the data on the destination grid should be written

        varlist : list of str, optional
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

        Raises
        ------
        OSError
            If ``ncremap`` is not in the system path.

        ValueError
            If ``mapping_filename`` is ``None`` (meaning no remapping is
            needed).
        """

        if self.mapping_filename is None:
            raise ValueError('No mapping file was given because remapping is '
                             'not necessary. The calling\n'
                             'code should simply use the constents of {} '
                             'directly.'.format(infilename))

        if not overwrite and os.path.exists(outfilename):
            # a remapped file already exists, so nothing to do
            return

        if isinstance(self.src_descrip, PointCollectionDescriptor):
            raise TypeError('Source grid is a point collection, which is not'
                            'supported.')

        if find_executable('ncremap') is None:
            raise OSError('ncremap not found. Make sure the latest nco '
                          'package is installed: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        tempfiles = []

        infilename, unpermute_dims = self._permute_dimensions(infilename)
        if bool(unpermute_dims):
            # we permuted som dimensions and saved to a temp file
            tempfiles.append(infilename)

        regird_args = []

        src_coords = self.src_descrip.lon_lat_coords
        if src_coords is not None:
            regird_args.extend(
                ['--rgr lat_nm={}'.format(src_coords[0]),
                 '--rgr lon_nm={}'.format(src_coords[1])])

        dst_coords = self.dst_descrip.lon_lat_coords
        if dst_coords is not None:
            regird_args.extend(
                ['--rgr lat_nm_out={}'.format(dst_coords[0]),
                 '--rgr lon_nm_out={}'.format(dst_coords[1])])

        dst_dims = list(self.dst_descrip.sizes.keys())
        if len(dst_dims) == 2:
            regird_args.extend(
                ['--rgr lat_dmn_nm={}'.format(dst_dims[0]),
                 '--rgr lon_dmn_nm={}'.format(dst_dims[1])])
        elif len(dst_dims) == 1:
            regird_args.extend(
                ['--rgr col_nm={}'.format(dst_dims[0])])

        add_coords = False
        if self.dst_descrip.coords is not None:
            if src_coords is None:
                add_coords = True
            else:
                add_coords = any([name not in src_coords for name in
                                  self.dst_descrip.coords.keys()])

        if add_coords or bool(unpermute_dims):
            remap_filename = _get_temp_path()
            tempfiles.append(remap_filename)
        else:
            remap_filename = outfilename

        args = ['ncremap',
                '-i', infilename,
                '-m', self.mapping_filename,
                '--vrb=1',
                '-o', remap_filename]

        if renormalize is not None:
            args.append('--renormalization_threshold={}'.format(renormalize))

        if len(regird_args) > 0:
            args.extend(['-R', ' '.join(regird_args)])

        if varlist is not None:
            args.extend(['-v', ','.join(varlist)])

        # set an environment variable to make sure we're not using czender's
        # local version of NCO instead of one we have intentionally loaded
        env = os.environ.copy()
        env['NCO_PATH_OVERRIDE'] = 'No'

        if logger is None:
            print('running: {}'.format(' '.join(args)))
            # make sure any output is flushed before we add output from the
            # subprocess
            sys.stdout.flush()
            sys.stderr.flush()

            subprocess.check_call(args, env=env)
        else:
            logger.info('running: {}'.format(' '.join(args)))
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

        if bool(unpermute_dims):
            # we permuted som dimensions, so we need to unpermute them
            if add_coords:
                unpermute_filename = _get_temp_path()
                tempfiles.append(unpermute_filename)
            else:
                unpermute_filename = outfilename
            self._unpermute_dimensions(remap_filename, unpermute_filename,
                                       unpermute_dims)
            remap_filename = unpermute_filename

        if add_coords:
            # add missing x and y coordinates
            ds = xarray.open_dataset(remap_filename)
            # load into memory because we're going to write to the same file
            ds.load()
            ds.close()

            for name, coord in self.dst_descrip.coords.items():
                # should overwrite coords if they already exist
                ds.coords[name] = xarray.DataArray.from_dict(coord)
            write_netcdf(ds, outfilename)

        for filename in tempfiles:
            os.remove(filename)

    def remap(self, ds, renorm_thresh=None):
        """
        Given a source data set, returns a remapped version of the data set,
        possibly masked and renormalized.

        Parameters
        ----------
        ds : ``xarray.Dataset`` or ``xarray.DataArray``
            The dimention(s) along ``self.sourceDimNames`` must match
            ``self.src_grid_dims`` read from the mapping file.

        renorm_thresh : float, optional
            The minimum weight of a denstination cell after remapping, below
            which it is masked out, or ``None`` for no renormalization and
            masking.

        Returns
        -------
        ds_remap : `xarray.Dataset`` or ``xarray.DataArray``
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

        if self.mapping_filename is None:
            # No remapping is needed
            return ds

        self._load_mapping()

        for dim, size in self.src_descrip.sizes.items():
            if size != ds.sizes[dim]:
                raise ValueError('data set and remapping source dimension {} '
                                 'don\'t have the same size: {} != {}'.format(
                                     dim, size,
                                     ds.sizes[dim]))

        if isinstance(ds, xarray.DataArray):
            ds_remap = self._remap_data_array(ds, renorm_thresh)
        elif isinstance(ds, xarray.Dataset):
            drop = []
            for var in ds.data_vars:
                if self._check_drop(ds[var]):
                    drop.append(var)
            ds_remap = ds.drop_vars(drop)
            ds_remap = ds_remap.map(self._remap_data_array,
                                    keep_attrs=True,
                                    args=(renorm_thresh,))

            for name, coord in self.dst_descrip.coords.items():
                if name not in ds_remap.coords:
                    ds_remap.coords[name] = xarray.DataArray.from_dict(coord)

        else:
            raise TypeError('ds not an xarray Dataset or DataArray.')

        # Update history attribute of netCDF file
        if 'history' in ds_remap.attrs:
            newhist = '\n'.join([ds_remap.attrs['history'],
                                 ' '.join(sys.argv[:])])
        else:
            newhist = ' '.join(sys.argv[:])
        ds_remap.attrs['history'] = newhist

        ds_remap.attrs['meshname'] = self.dst_descrip.meshname

        return ds_remap

    def _load_mapping(self):
        """
        Load weights and indices from a mapping file, if this has not already
        been done
        """

        if self._mapping_loaded:
            return

        ds_mapping = xarray.open_dataset(self.mapping_filename)
        n_a = ds_mapping.dims['n_a']
        n_b = ds_mapping.dims['n_b']

        src_dim_count = len(self.src_descrip.sizes)
        src_grid_rank = ds_mapping.dims['src_grid_rank']
        dst_dim_count = len(self.dst_descrip.sizes)
        dst_grid_rank = ds_mapping.dims['dst_grid_rank']

        # check that the mapping file has the right number of dimensions
        if src_dim_count != src_grid_rank or \
                dst_dim_count != dst_grid_rank:
            raise ValueError('The number of source and/or '
                             'destination dimensions does not\n'
                             'match the expected number of source and '
                             'destination dimensions in the mapping\n'
                             'file. {} != {} and/or {} != {}'.format(
                                 src_dim_count, src_grid_rank,
                                 dst_dim_count, dst_grid_rank))

        # grid dimensions need to be reversed because they are in Fortran order
        self.src_grid_dims = ds_mapping['src_grid_dims'].values[::-1]
        self.dst_grid_dims = ds_mapping['dst_grid_dims'].values[::-1]

        # now, check that each source and destination dimension is right
        Remapper._check_dims(self.src_descrip.sizes,  self.src_grid_dims)
        Remapper._check_dims(self.dst_descrip.sizes,  self.dst_grid_dims)

        self.frac_b = ds_mapping['frac_b'].values

        col = ds_mapping['col'].values - 1
        row = ds_mapping['row'].values - 1
        S = ds_mapping['S'].values
        self.matrix = csr_matrix((S, (row, col)), shape=(n_b, n_a))

        self._mapping_loaded = True

    def _permute_dimensions(self, infilename):
        src_dims = list(self.src_descrip.sizes.keys())
        dst_dims = list(self.dst_descrip.sizes.keys())
        ds = xarray.open_dataset(infilename)
        varnames = list(ds.data_vars.keys()) + list(ds.coords.keys())

        # skip variables or coordinates that are going to be replaced by output
        # coordinates
        varnames = [name for name in varnames if name not in
                    self.dst_descrip.coords]

        unpermute_dims = {}
        for name in varnames:
            var = ds[name]
            orig = list(var.dims)
            new = list(orig)
            unpermute = list(orig)
            dims_in_var = all([dim in orig for dim in src_dims])

            if dims_in_var:
                addIndex = len(orig)
                for dim in src_dims:
                    if dim in new:
                        index = new.index(dim)
                        # move it to the end
                        new.append(new.pop(index))
                        unpermute.pop(index)
                        addIndex = min(addIndex, index)

                if new != orig:
                    # we've permuted dims
                    unpermute[addIndex:addIndex] = dst_dims
                    unpermute_dims[name] = unpermute
                    ds[name] = var.transpose(*new, transpose_coords=True)

        if bool(unpermute_dims):
            # we add something to the unpermuted dimensions dict, so we need to
            # write out the permuted
            permuted_filename = _get_temp_path()
            write_netcdf(ds, permuted_filename)
        else:
            permuted_filename = infilename

        return permuted_filename, unpermute_dims

    # noinspection PyMethodMayBeStatic
    def _unpermute_dimensions(self, remap_filename, unpermute_filename,
                              unpermute_dims):
        ds = xarray.open_dataset(remap_filename)
        for name, dims in unpermute_dims.items():
            if name in ds:
                ds[name] = ds[name].transpose(*dims, transpose_coords=True)
        write_netcdf(ds, unpermute_filename)

    @staticmethod
    def _check_dims(sizes, grid_dims):
        for index, (dim, dim_size) in enumerate(sizes.items()):
            check_dim_size = grid_dims[index]
            if dim_size != check_dim_size:
                raise ValueError('mesh descriptor and remapping '
                                 'dimension {} don\'t have the same size: \n'
                                 '{} != {}'.format(dim, dim_size,
                                                   check_dim_size))

    def _check_drop(self, da):
        src_dims = self.src_descrip.sizes.keys()

        src_dims_in_array = [dim in da.dims for dim in src_dims]

        return (numpy.any(src_dims_in_array) and not
                numpy.all(src_dims_in_array))

    def _remap_data_array(self, da, renorm_thresh):
        """
        Remap a single xarray data array
        """

        src_dims = self.src_descrip.sizes.keys()
        dst_dims = self.dst_descrip.sizes.keys()

        src_dims_in_array = [dim in da.dims for dim in src_dims]

        if not numpy.any(src_dims_in_array):
            # no remapping is needed
            return da

        if not numpy.all(src_dims_in_array):
            # no remapping is possible so the variable array should have been
            # dropped
            raise ValueError('Data array with some (but not all) required '
                             'source dims cannot be remapped\n'
                             'and should have been dropped.')

        # make a list of dims and remap_axes
        dims = []
        remap_axes = []
        dst_dims_added = False
        for index, dim in enumerate(da.dims):
            if dim in src_dims:
                remap_axes.append(index)
                if not dst_dims_added:
                    dims.extend(dst_dims)
                    dst_dims_added = True
            else:
                dims.append(dim)

        # make a dict of coords
        coord_dict = {}
        # copy unmodified coords
        for coord in da.coords:
            src_dim_in_coord = numpy.any([dim in da.coords[coord].dims
                                          for dim in src_dims])
            if not src_dim_in_coord:
                coord_dict[coord] = {'dims': da.coords[coord].dims,
                                     'data': da.coords[coord].values,
                                     'attrs': da.coords[coord].attrs}

        # add dest coords
        for name, coord in self.dst_descrip.coords.items():
            coord_dims = coord['dims']
            if all([dim in dims for dim in coord_dims]):
                coord_dict[name] = coord

        # remap the values
        field = da.values
        mask = numpy.isnan(field)
        if numpy.count_nonzero(mask) > 0:
            field = numpy.ma.masked_array(field, mask)
        remapped_field = self._remap_numpy_array(field, remap_axes,
                                                 renorm_thresh)

        array_dict = {'coords': coord_dict,
                      'attrs': da.attrs,
                      'dims': dims,
                      'data': remapped_field,
                      'name': da.name}

        # make a new data array
        da_remap = xarray.DataArray.from_dict(array_dict)

        return da_remap

    def _remap_numpy_array(self, in_field, remap_axes,
                           renorm_thresh):
        """
        Remap a single numpy array
        """

        # permute the dimensions of in_field so the axes to remap are first,
        # then flatten the remapping and the extra dimensions separately for
        # the matrix multiply
        extra_axes = [axis for axis in numpy.arange(in_field.ndim)
                      if axis not in remap_axes]

        new_shape = [numpy.prod([in_field.shape[axis] for axis in remap_axes])]
        if len(extra_axes) > 0:
            extra_shape = [in_field.shape[axis] for axis in extra_axes]
            new_shape.append(numpy.prod(extra_shape))
        else:
            extra_shape = []
            new_shape.append(1)

        permuted_axes = remap_axes + extra_axes

        # permute axes so the remapped dimension(s) come first and "flatten"
        # the remapping dimension
        in_field = in_field.transpose(permuted_axes).reshape(new_shape)

        masked = isinstance(in_field, numpy.ma.MaskedArray)
        threshold = renorm_thresh is not None
        if threshold:
            if masked:
                in_mask = numpy.array(numpy.logical_not(in_field.mask), float)
                out_field = self.matrix.dot(in_mask * in_field)
                out_mask = self.matrix.dot(in_mask)
                mask = out_mask > renorm_thresh
            else:
                out_field = self.matrix.dot(in_field)
                # make frac_b match the shape of out_field
                out_mask = numpy.reshape(self.frac_b,
                                         (len(self.frac_b), 1)).repeat(
                    new_shape[1], axis=1)
                mask = out_mask > renorm_thresh

            # normalize the result based on out_mask
            out_field[mask] /= out_mask[mask]
            out_field = numpy.ma.masked_array(out_field,
                                              mask=numpy.logical_not(mask))
        else:
            if masked:
                in_mask = numpy.array(numpy.logical_not(in_field.mask), float)
                out_field = self.matrix.dot(in_mask * in_field)
            else:
                out_field = self.matrix.dot(in_field)

        dst_remap_dim_count = len(self.dst_grid_dims)
        out_dim_count = len(extra_shape) + dst_remap_dim_count

        # "unflatten" the remapped dimension(s)
        dst_shape = list(self.dst_grid_dims) + extra_shape
        out_field = numpy.reshape(out_field, dst_shape)

        # "unpermute" the axes to be in the expected order
        index = numpy.amin(remap_axes)
        unpermute_axes = list(numpy.arange(dst_remap_dim_count, out_dim_count))
        unpermute_axes = (unpermute_axes[0:index] +
                          list(numpy.arange(dst_remap_dim_count)) +
                          unpermute_axes[index:])
        out_field = numpy.transpose(out_field, axes=unpermute_axes)

        return out_field


# noinspection PyProtectedMember
def _get_temp_path():
    """Returns the name of a temporary NetCDF file"""
    return '{}/{}.nc'.format(tempfile._get_default_tempdir(),
                             next(tempfile._get_candidate_names()))
