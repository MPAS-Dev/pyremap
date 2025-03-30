import os
from tempfile import TemporaryDirectory

from pyremap.remapper.setup import _setup_remapper
from pyremap.utility import check_call


def _build_map(remapper):
    """
    Make the mapping file
    """
    _setup_remapper(remapper)
    map_tool = remapper.map_tool

    src_scrip_filename = remapper.src_scrip_filename
    dst_scrip_filename = remapper.dst_scrip_filename
    scrip_dir = os.path.dirname(src_scrip_filename)

    if remapper.use_tmp:
        tempobj = TemporaryDirectory()
        scrip_dir = os.path.join(tempobj.name, scrip_dir)
        src_scrip_filename = os.path.join(
            tempobj.name, src_scrip_filename
        )
        dst_scrip_filename = os.path.join(
            tempobj.name, dst_scrip_filename
        )
    else:
        tempobj = None

    if scrip_dir != '':
        os.makedirs(scrip_dir, exist_ok=True)

    src_descriptor = remapper.src_descriptor
    src_descriptor.to_scrip(src_scrip_filename)

    dst_descriptor = remapper.dst_descriptor
    dst_descriptor.to_scrip(
        dst_scrip_filename,
        expand_dist=remapper.expand_dist,
        expand_factor=remapper.expand_factor,
    )

    if map_tool == 'esmf':
        args = _esmf_build_map_args(
            remapper, src_scrip_filename, dst_scrip_filename)

        esmf_path = remapper.esmf_path
        if esmf_path is None:
            esmf_exe = 'ESMF_RegridWeightGen'
        else:
            esmf_exe = os.path.join(
                esmf_path, 'bin', 'ESMF_RegridWeightGen')
        args = [esmf_exe] + args

    elif map_tool == 'moab':
        moab_path = remapper.moab_path
        src_scrip_filename = _moab_partition_scrip_file(
            remapper, src_scrip_filename, moab_path
        )
        dst_scrip_filename = _moab_partition_scrip_file(
            remapper, dst_scrip_filename, moab_path
        )
        args = _moab_build_map_args(
            remapper, src_scrip_filename, dst_scrip_filename)

        if moab_path is None:
            moab_exe = 'mbtempest'
        else:
            moab_exe = os.path.join(moab_path, 'bin', 'mbtempest')
        args = [moab_exe] + args

    ntasks = remapper.ntasks
    if ntasks > 1:
        parallel_exec = remapper.parallel_exec
        if parallel_exec == 'srun':
            args = [parallel_exec, '-n', str(ntasks)] + args
        elif parallel_exec == 'mpirun':
            args = [parallel_exec, '-np', str(ntasks)] + args
        else:
            raise ValueError(
                f'Unexpected parallel_exec {parallel_exec}. '
                f'Valid values are "mpirun" or "srun".'
            )

    check_call(args)

    if tempobj is not None:
        tempobj.cleanup()


def _moab_partition_scrip_file(remapper, in_filename, moab_path):
    """
    Partition SCRIP file for parallel mbtempest use
    """
    ntasks = remapper.ntasks

    print(f'Partition SCRIP file {in_filename}')

    h5m_filename = in_filename.replace('.nc', '.h5m')
    h5m_part_filename = in_filename.replace('.nc', f'.p{ntasks}.h5m')

    if moab_path is None:
        mbconvert = 'mbconvert'
        mbpart = 'mbpart'
    else:
        mbconvert = os.path.join(moab_path, 'bin', 'mbconvert')
        mbpart = os.path.join(moab_path, 'bin', 'mbpart')

    # Convert source SCRIP to mbtempest
    args = [
        mbconvert,
        '-B',
        in_filename,
        h5m_filename,
    ]
    check_call(args)

    # Partition source SCRIP
    args = [
        mbpart,
        f'{ntasks}',
        '-z',
        'RCB',
        h5m_filename,
        h5m_part_filename,
    ]
    check_call(args)

    print('  Done.')

    return h5m_part_filename


def _esmf_build_map_args(remapper, src_scrip_filename, dst_scrip_filename):
    """
    Get command-line arguments for making a mapping file with
    ESMF_RegridWeightGen
    """

    args = [
        '--source',
        src_scrip_filename,
        '--destination',
        dst_scrip_filename,
        '--weight',
        remapper.map_filename,
        '--method',
        remapper.method,
        '--netcdf4',
    ]

    if remapper.src_descriptor.regional:
        args.append('--src_regional')

    if remapper.dst_descriptor.regional:
        args.append('--dst_regional')

    if remapper.src_descriptor.regional or remapper.dst_descriptor.regional:
        args.append('--ignore_unmapped')

    return args


def _moab_build_map_args(remapper, src_scrip_filename, dst_scrip_filename):
    """
    Get command-line arguments for making a mapping file with mbtempest
    """
    fvmethod = {'conserve': 'none', 'bilinear': 'bilin'}

    args = [
        '--type',
        '5',
        '--load',
        src_scrip_filename,
        '--load',
        dst_scrip_filename,
        '--file',
        remapper.map_filename,
        '--weights',
        '--gnomonic',
        '--boxeps',
        '1e-9',
        '--method',
        'fv',
        '--method',
        'fv',
        '--order',
        '1',
        '--order',
        '1',
        '--fvmethod',
        fvmethod[remapper.method],
    ]

    return args
