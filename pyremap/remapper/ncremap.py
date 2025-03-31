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
    """
    Given a source file defining either an MPAS mesh or a lat-lon grid and
    a destination file or set of arrays defining a lat-lon grid, constructs
    a mapping file used for interpolation between the source and
    destination grids."
    """
    map_filename = remapper.map_filename
    src_descriptor = remapper.src_descriptor
    dst_descriptor = remapper.dst_descriptor
    if map_filename is None:
        raise ValueError('No mapping file has been defined')

    if not overwrite and os.path.exists(out_filename):
        # a remapped file already exists, so nothing to do
        return

    if isinstance(src_descriptor, PointCollectionDescriptor):
        raise TypeError(
            'Source grid is a point collection, which is notsupported.'
        )

    if which('ncremap') is None:
        raise OSError(
            'ncremap not found. Make sure NCO is available in your PATH.'
        )

    if parallel_exec is not None:
        # use the specified parallel executable
        args = parallel_exec.split(' ')
    else:
        args = []

    args.extend(['ncremap', '-m', map_filename, '--vrb=1'])

    regrid_args = []

    mpas_descriptors = (
        MpasCellMeshDescriptor,
        MpasEdgeMeshDescriptor,
        MpasVertexMeshDescriptor,
    )
    if isinstance(src_descriptor, mpas_descriptors):
        args.extend(['-P', 'mpas'])
        if not replace_mpas_fill:
            # the -C (climatology) flag prevents ncremap from trying to
            # add a _FillValue attribute that might already be present
            # and quits with an error
            args.append('-C')

    if isinstance(src_descriptor, MpasEdgeMeshDescriptor):
        regrid_args.extend(['--rgr col_nm=nEdges'])

    if (
        isinstance(src_descriptor, mpas_descriptors)
        and renormalize is not None
    ):
        # we also want to make sure cells that receive no data are
        # marked with fill values, even if the source MPAS data
        # doesn't have a fill value
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
    if isinstance(dst_descriptor, PointCollectionDescriptor):
        regrid_args.extend(['--rgr lat_nm_out=lat', '--rgr lon_nm_out=lon'])

    if len(regrid_args) > 0:
        args.extend(['-R', ' '.join(regrid_args)])

    # set an environment variable to make sure we're not using czender's
    # local version of NCO instead of one we have intentionally loaded
    env = os.environ.copy()
    env['NCO_PATH_OVERRIDE'] = 'No'

    args.extend([in_filename, out_filename])

    check_call(args, logger=logger, env=env)
