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

import sys

import numpy as np
import xarray as xr
from scipy.sparse import csr_matrix


def _remap_numpy(remapper, ds, renormalization_threshold):
    """
    Given a source data set, returns a remapped version of the data set,
    possibly masked and renormalized.
    """
    map_filename = remapper.map_filename
    if map_filename is None:
        raise ValueError('No mapping file has been defined')

    _load_mapping(remapper)

    src_grid_dims = remapper._ds_map['src_grid_dims'].values[::-1]

    for index, dim in enumerate(remapper.src_descriptor.dims):
        if src_grid_dims[index] != ds.sizes[dim]:
            raise ValueError(
                f"data set and remapping source dimension {dim} don't "
                f'have the same size: {src_grid_dims[index]} != '
                f'{ds.sizes[dim]}'
            )

    if isinstance(ds, xr.DataArray):
        ds_remap = _remap_data_array(ds, remapper, renormalization_threshold)
    elif isinstance(ds, xr.Dataset):
        drop = []
        for var in ds.data_vars:
            if _check_drop(remapper, ds[var]):
                drop.append(var)
        ds_remap = ds.drop_vars(drop)
        ds_remap = ds_remap.map(
            _remap_data_array,
            keep_attrs=True,
            args=(
                remapper,
                renormalization_threshold,
            ),
        )
    else:
        raise TypeError('ds not an xarray Dataset or DataArray.')

    # Update history attribute of netCDF file
    current_hist = ' '.join(sys.argv[:])
    if 'history' in ds_remap.attrs:
        newhist = '\n'.join([ds_remap.attrs['history'], current_hist])
    else:
        newhist = current_hist
    ds_remap.attrs['history'] = newhist

    ds_remap.attrs['mesh_name'] = remapper.dst_descriptor.mesh_name

    return ds_remap


def _load_mapping(remapper):
    """
    Load weights and indices from a mapping file, if this has not already
    been done
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if remapper._ds_map is not None:
        return

    src_descriptor = remapper.src_descriptor
    dst_descriptor = remapper.dst_descriptor
    map_filename = remapper.map_filename

    ds_map = xr.open_dataset(map_filename)
    n_a = ds_map.sizes['n_a']
    n_b = ds_map.sizes['n_b']

    n_source_dims = len(src_descriptor.dims)
    src_grid_rank = ds_map.sizes['src_grid_rank']
    n_destination_dims = len(dst_descriptor.dims)
    dst_grid_rank = ds_map.sizes['dst_grid_rank']

    # check that the mapping file has the right number of dimensions
    if n_source_dims != src_grid_rank or n_destination_dims != dst_grid_rank:
        raise ValueError(
            f'The number of source and/or destination dimensions does not '
            f'match the expected \n'
            f'number of source and destination dimensions in the mapping '
            f'file. \n'
            f'{n_source_dims} != {src_grid_rank} and/or {n_destination_dims} '
            f'!= {dst_grid_rank}'
        )

    # grid dimensions need to be reversed because they are in Fortran order
    src_grid_dims = ds_map['src_grid_dims'].values[::-1]
    dst_grid_dims = ds_map['dst_grid_dims'].values[::-1]

    # now, check that each source and destination dimension is right
    for index in range(len(src_descriptor.dims)):
        dim = src_descriptor.dims[index]
        dim_size = src_descriptor.dim_sizes[index]
        check_dim_size = src_grid_dims[index]
        if dim_size != check_dim_size:
            raise ValueError(
                f'source mesh descriptor and remapping source dimension '
                f"{dim} don't have the same size: \n"
                f'{dim_size} != {check_dim_size}'
            )
    for index in range(len(dst_descriptor.dims)):
        dim = dst_descriptor.dims[index]
        dim_size = dst_descriptor.dim_sizes[index]
        check_dim_size = dst_grid_dims[index]
        if dim_size != check_dim_size:
            raise ValueError(
                f'dest. mesh descriptor and remapping dest. dimension '
                f"{dim} don't have the same size: \n"
                f'{dim_size} != {check_dim_size}'
            )

    col = ds_map['col'].values - 1
    row = ds_map['row'].values - 1
    s = ds_map['S'].values
    remapper._matrix = csr_matrix((s, (row, col)), shape=(n_b, n_a))

    remapper._ds_map = ds_map


def _check_drop(remapper, da):
    src_dims = remapper.src_descriptor.dims

    src_dims_in_array = [dim in da.dims for dim in src_dims]

    return np.any(src_dims_in_array) and not np.all(src_dims_in_array)


def _remap_data_array(da, remapper, renormalization_threshold):
    """
    Remap a single xarray data array
    """
    src_dims = remapper.src_descriptor.dims
    dst_dims = remapper.dst_descriptor.dims

    src_dims_in_array = [dim in da.dims for dim in src_dims]

    if not np.any(src_dims_in_array):
        # no remapping is needed
        return da

    if not np.all(src_dims_in_array):
        # no remapping is possible so the variable array should have been
        # dropped
        raise ValueError(
            'Data array with some (but not all) required source dims cannot '
            'be remapped and should have been dropped.'
        )

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
        src_dim_in_coord = np.any(
            [dim in da.coords[coord].dims for dim in src_dims]
        )
        if not src_dim_in_coord:
            coord_dict[coord] = {
                'dims': da.coords[coord].dims,
                'data': da.coords[coord].values,
            }

    # add destination coords
    coord_dict.update(remapper.dst_descriptor.coords)

    # remap the values
    field = da.values
    mask = np.isnan(field)
    if np.count_nonzero(mask) > 0:
        field = np.ma.masked_array(field, mask)
    remapped_field = _remap_numpy_array(
        remapper, field, remap_axes, renormalization_threshold
    )

    array_dict = {
        'coords': coord_dict,
        'attrs': da.attrs,
        'dims': dims,
        'data': remapped_field,
        'name': da.name,
    }

    # make a new data array
    remapped_array = xr.DataArray.from_dict(array_dict)

    return remapped_array


def _remap_numpy_array(
    remapper, in_field, remap_axes, renormalization_threshold
):
    """
    Remap a single numpy array
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    # permute the dimensions of in_field so the axes to remap are first,
    # then flatten the remapping and the extra dimensions separately for
    # the matrix multiply
    extra_axes = [
        axis for axis in np.arange(in_field.ndim) if axis not in remap_axes
    ]

    new_shape = [np.prod([in_field.shape[axis] for axis in remap_axes])]
    if len(extra_axes) > 0:
        extra_shape = [in_field.shape[axis] for axis in extra_axes]
        new_shape.append(np.prod(extra_shape))
    else:
        extra_shape = []
        new_shape.append(1)

    permuted_axes = remap_axes + extra_axes

    matrix = remapper._matrix
    # grid dimensions need to be reversed because they are in Fortran order
    dst_grid_dims = remapper._ds_map['dst_grid_dims'].values[::-1]

    # permute axes so the remapped dimension(s) come first and "flatten"
    # the remapping dimension
    in_field = in_field.transpose(permuted_axes).reshape(new_shape)

    masked = (
        isinstance(in_field, np.ma.MaskedArray)
        and renormalization_threshold is not None
    )
    if masked:
        in_mask = np.array(np.logical_not(in_field.mask), float)
        out_field = matrix.dot(in_mask * in_field)
        out_mask = matrix.dot(in_mask)
        mask = out_mask > renormalization_threshold
    else:
        out_field = matrix.dot(in_field)
        # make frac_b match the shape of out_field
        frac_b = remapper._ds_map['frac_b'].values
        out_mask = np.reshape(frac_b, (len(frac_b), 1)).repeat(
            new_shape[1], axis=1
        )
        mask = out_mask > 0.0

    # normalize the result based on out_mask
    out_field[mask] /= out_mask[mask]
    out_field = np.ma.masked_array(out_field, mask=np.logical_not(mask))

    dest_remap_dim_count = len(dst_grid_dims)
    out_dim_count = len(extra_shape) + dest_remap_dim_count

    # "unflatten" the remapped dimension(s)
    dest_shape = list(dst_grid_dims) + extra_shape
    out_field = np.reshape(out_field, dest_shape)

    # "unpermute" the axes to be in the expected order
    index = np.amin(remap_axes)
    unpermute_axes = list(np.arange(dest_remap_dim_count, out_dim_count))
    unpermute_axes = (
        unpermute_axes[0:index]
        + list(np.arange(dest_remap_dim_count))
        + unpermute_axes[index:]
    )
    out_field = np.transpose(out_field, axes=unpermute_axes)

    return out_field
