import os
from tempfile import TemporaryDirectory

from pyremap.remapper.setup import _setup_remapper
from pyremap.utility import check_call


def _build_map(remapper, logger=None):
    """
    Make the mapping file

    Parameters
    ----------
    remapper : pyremap.remapper.MpasRemapper
        The remapper object to use to create the mapping file
    logger : logging.Logger, optional
        A logger to which ncclimo output should be redirected
    """
    if remapper.use_tmp:
        with TemporaryDirectory() as tmp_dir:
            args = _write_files_and_get_args(remapper, tmp_dir, logger)
            _run_map_tool(remapper, args, logger)
    else:
        args = _write_files_and_get_args(remapper, '', logger)
        _run_map_tool(remapper, args, logger)


def _prepare_map(remapper, logger=None):
    """
    Write the files the mapping tool reads and get the command that makes
    the mapping file, without running it

    Parameters
    ----------
    remapper : pyremap.remapper.MpasRemapper
        The remapper object to use to create the mapping file
    logger : logging.Logger, optional
        A logger to which output should be redirected

    Returns
    -------
    args : list of str
        The mapping tool and its arguments, with no parallel launcher
    """
    if remapper.use_tmp:
        raise ValueError(
            'prepare_map() requires use_tmp=False, since the mapping tool '
            'reads the files it writes after it returns.'
        )
    return _write_files_and_get_args(remapper, '', logger)


def _write_files_and_get_args(remapper, tmp_dir, logger):
    """
    Write the SCRIP files (and, for MOAB, the HDF5 files made from them) in
    ``tmp_dir`` and get the command that makes the mapping file
    """
    _setup_remapper(remapper)
    map_tool = remapper.map_tool

    src_scrip_filename = os.path.join(tmp_dir, remapper.src_scrip_filename)
    dst_scrip_filename = os.path.join(tmp_dir, remapper.dst_scrip_filename)
    scrip_dir = os.path.dirname(src_scrip_filename)

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
            remapper, src_scrip_filename, dst_scrip_filename
        )

        esmf_path = remapper.esmf_path
        if esmf_path is None:
            esmf_exe = 'ESMF_RegridWeightGen'
        else:
            esmf_exe = os.path.join(esmf_path, 'bin', 'ESMF_RegridWeightGen')
        args = [esmf_exe] + args

    elif map_tool == 'moab':
        moab_path = remapper.moab_path
        src_scrip_filename = _moab_convert_scrip_file(
            remapper, src_scrip_filename, moab_path, logger
        )
        dst_scrip_filename = _moab_convert_scrip_file(
            remapper, dst_scrip_filename, moab_path, logger
        )
        args = _moab_build_map_args(
            remapper, src_scrip_filename, dst_scrip_filename
        )

        if moab_path is None:
            moab_exe = 'mbtempest'
        else:
            moab_exe = os.path.join(moab_path, 'bin', 'mbtempest')
        args = [moab_exe] + args

    return args


def _run_map_tool(remapper, args, logger):
    """
    Run the mapping tool, with the parallel launcher if there is more than
    one task
    """
    ntasks = remapper.ntasks
    if ntasks > 1:
        parallel_exec = remapper.parallel_exec
        parallel_args = parallel_exec.split(' ')
        if parallel_args[0] == 'srun':
            args = parallel_args + ['-n', str(ntasks)] + args
        elif parallel_args[0] == 'mpirun':
            args = parallel_args + ['-np', str(ntasks)] + args
        else:
            raise ValueError(
                f'Unexpected parallel_exec {parallel_exec}. '
                f'Valid values are "mpirun" or "srun" with optional flags.'
            )

    check_call(args, logger=logger)


def _moab_convert_scrip_file(remapper, in_filename, moab_path, logger):
    """
    Convert SCRIP file to HDF5 for mbtempest use, partitioning it if more than
    one MPI task will be used
    """
    ntasks = remapper.ntasks

    print(f'Convert SCRIP file {in_filename}')

    h5m_filename = in_filename.replace('.nc', '.h5m')

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
    check_call(args, logger=logger)

    if ntasks == 1:
        # nothing to partition (mbpart requires more than one part), so
        # mbtempest reads the unpartitioned mesh.  The conversion above is
        # still needed: mbtempest fails to compute weights when a mesh with
        # more than 4 vertices per cell (e.g. an MPAS cell mesh) is loaded
        # directly from SCRIP as the destination mesh.
        print('  Done.')
        return h5m_filename

    h5m_part_filename = in_filename.replace('.nc', f'.p{ntasks}.h5m')

    # Partition source SCRIP
    args = [
        mbpart,
        f'{ntasks}',
        '-z',
        'RCB',
        h5m_filename,
        h5m_part_filename,
    ]
    check_call(args, logger=logger)

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
