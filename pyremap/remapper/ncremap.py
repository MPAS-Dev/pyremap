import os
from shutil import which

from pyremap.descriptor import (
    LatLonGridDescriptor,
    MpasCellMeshDescriptor,
    MpasEdgeMeshDescriptor,
    MpasVertexMeshDescriptor,
    PointCollectionDescriptor,
    ProjectionGridDescriptor,
)
from pyremap.utility import check_call


def _validate_inputs(remapper, out_filename, overwrite):
    if remapper.map_filename is None:
        raise ValueError('No mapping file has been defined')
    if not overwrite and os.path.exists(out_filename):
        return False  # Skip processing
    if isinstance(remapper.src_descriptor, PointCollectionDescriptor):
        raise TypeError(
            'Source grid is a point collection, which is not supported.'
        )
    if which('ncremap') is None:
        raise OSError(
            'ncremap not found. Make sure NCO is available in your PATH.'
        )
    return True


def _build_args(
    remapper, variable_list, replace_mpas_fill, renormalize, parallel_exec
):
    args = parallel_exec.split(' ') if parallel_exec else []
    args.extend(['ncremap', '-m', remapper.map_filename, '--vrb=1'])

    src_descriptor = remapper.src_descriptor
    dst_descriptor = remapper.dst_descriptor
    regrid_args = []

    mpas_descriptors = (
        MpasCellMeshDescriptor,
        MpasEdgeMeshDescriptor,
        MpasVertexMeshDescriptor,
    )
    oned_descriptors = (
        MpasCellMeshDescriptor,
        MpasEdgeMeshDescriptor,
        MpasVertexMeshDescriptor,
        PointCollectionDescriptor,
    )

    if isinstance(src_descriptor, mpas_descriptors):
        args.extend(['-P', 'mpas'])
        if not replace_mpas_fill:
            args.append('-C')
    if isinstance(src_descriptor, MpasEdgeMeshDescriptor):
        regrid_args.extend(['--rgr col_nm=nEdges'])
    if (
        isinstance(src_descriptor, mpas_descriptors)
        and renormalize is not None
    ):
        args.append('--add_fill_value')
    if variable_list is not None:
        args.extend(['-v', ','.join(variable_list)])
    if renormalize is not None:
        regrid_args.append(f'--renormalize={renormalize}')
    if isinstance(src_descriptor, LatLonGridDescriptor):
        regrid_args.extend(
            [
                f'--rgr lat_nm={src_descriptor.lat_var_name}',
                f'--rgr lon_nm={src_descriptor.lon_var_name}',
            ]
        )
    elif isinstance(src_descriptor, ProjectionGridDescriptor):
        regrid_args.extend(
            [
                f'--rgr lat_nm={src_descriptor.y_var_name}',
                f'--rgr lon_nm={src_descriptor.x_var_name}',
            ]
        )
    if isinstance(dst_descriptor, LatLonGridDescriptor):
        regrid_args.extend(
            [
                f'--rgr lat_nm_out={dst_descriptor.lat_var_name}',
                f'--rgr lon_nm_out={dst_descriptor.lon_var_name}',
                f'--rgr lat_dmn_nm={dst_descriptor.lat_var_name}',
                f'--rgr lon_dmn_nm={dst_descriptor.lon_var_name}',
            ]
        )
    elif isinstance(dst_descriptor, ProjectionGridDescriptor):
        regrid_args.extend(
            [
                f'--rgr lat_dmn_nm={dst_descriptor.y_var_name}',
                f'--rgr lon_dmn_nm={dst_descriptor.x_var_name}',
                '--rgr lat_nm_out=lat',
                '--rgr lon_nm_out=lon',
            ]
        )
    elif isinstance(dst_descriptor, oned_descriptors):
        dims = dst_descriptor.dims
        assert dims is not None, "dst descriptor's dims should not be None"
        regrid_args.extend([f'--rgr col_nm={dims[0]}'])
    if isinstance(dst_descriptor, PointCollectionDescriptor):
        regrid_args.extend(['--rgr lat_nm_out=lat', '--rgr lon_nm_out=lon'])
    if len(regrid_args) > 0:
        args.extend(['-R', ' '.join(regrid_args)])
    return args


def _set_environment():
    env = os.environ.copy()
    env['NCO_PATH_OVERRIDE'] = 'No'
    return env


def _ncremap(
    remapper,
    in_filename,
    out_filename,
    variable_list,
    overwrite,
    renormalize,
    logger,
    replace_mpas_fill,
    parallel_exec,
):
    # Validate inputs
    if not _validate_inputs(remapper, out_filename, overwrite):
        # file exists and overwrite is False
        return

    # Build arguments
    args = _build_args(
        remapper, variable_list, replace_mpas_fill, renormalize, parallel_exec
    )

    # Set environment variables
    env = _set_environment()

    # Add input and output filenames
    args.extend([in_filename, out_filename])

    # Execute the command
    check_call(args, logger=logger, env=env)
