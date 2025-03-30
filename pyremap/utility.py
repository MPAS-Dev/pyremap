import subprocess
import sys

import netCDF4
import numpy as np


def write_netcdf(
        ds, filename, format, engine=None, logger=None, fillvalues=None):
    """
    Write an xarray.Dataset to a file with NetCDF4 fill values.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to save

    filename : str
        The path for the NetCDF file to write

    format : {'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_DATA'}
        The NetCDF file format to use.

    engine : {'netcdf4', 'scipy', 'h5netcdf'}, optional
        The library to use for NetCDF output.  The default is the same as
        in :py:meth:`xarray.Dataset.to_netcdf` and depends on ``format``.

    fillvalues : dict, optional
        A dictionary of fill values for different NetCDF types.  Default is
        ``mpas_tools.io.default_fills``, which can be modified but which
        defaults to ``netCDF4.default_fillvals``
    """  # noqa: E501

    if fillvalues is None:
        fillvalues = netCDF4.default_fillvals

    encoding_dict = {}
    variable_names = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variable_name in variable_names:
        is_numeric = np.issubdtype(ds[variable_name].dtype, np.number)
        if is_numeric and np.any(np.isnan(ds[variable_name])):
            dtype = ds[variable_name].dtype
            for fill_type in fillvalues:
                if dtype == np.dtype(fill_type):
                    encoding_dict[variable_name] = \
                        {'_FillValue': fillvalues[fill_type]}
                    break
        else:
            encoding_dict[variable_name] = {'_FillValue': None}

    if format == 'NETCDF3_64BIT_DATA':
        # NETCDF3_64BIT_DATA needs special treatment because it isn't efficient
        # in xarray
        write_format = 'NETCDF4'
        write_filename = filename.replace('.nc', '_netcdf4.nc')
    else:
        write_format = format
        write_filename = filename

    ds.to_netcdf(
        write_filename,
        encoding=encoding_dict,
        format=write_format,
        engine=engine)

    if format == 'NETCDF3_64BIT_DATA':
        # Still need to convert to NETCDF3_64BIT_DATA
        args = ['ncks', '-O', '-5', write_filename, filename]
        check_call(args, logger=logger)


def check_call(args, logger=None, log_command=True, **kwargs):
    """
    Call the given subprocess

    Parameters
    ----------
    args : list or str
        A list or string of argument to the subprocess.  If ``args`` is a
        string, you must pass ``shell=True`` as one of the ``kwargs``.

    logger : logging.Logger, optional
        The logger to write output to

    log_command : bool, optional
        Whether to print the command that is running

    **kwargs : dict
        Keyword arguments to pass to subprocess.Popen

    Raises
    ------
    subprocess.CalledProcessError
        If the subprocess returns a non-zero exit code
    """

    if logger is None:
        # make sure any output is flushed before we add output from the
        # subprocess
        sys.stdout.flush()
        sys.stderr.flush()

        if log_command:
            _print_running(args, fn=print)
        subprocess.run(args, check=True, **kwargs)
    else:
        if log_command:
            _print_running(args, fn=logger.info)
        for handler in logger.handlers:
            handler.flush()

        process = subprocess.Popen(
            args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
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
            raise subprocess.CalledProcessError(
                process.returncode, ' '.join(args))


def _print_running(args, fn):
    if isinstance(args, str):
        print_args = args
    else:
        print_args = []
        for arg in args:
            if ' ' in arg:
                arg = f'"{arg}"'
            print_args.append(arg)
        print_args = ' '.join(print_args)
    fn(f'running: {print_args}')
