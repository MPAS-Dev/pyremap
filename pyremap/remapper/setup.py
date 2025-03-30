from pyremap import PointCollectionDescriptor
from pyremap.remapper.descriptor import _get_descriptor


def _setup_remapper(remapper):
    """
    Set up the descriptors and check the remapper
    """
    if remapper.src_descriptor is not None:
        src_descriptor = remapper.src_descriptor
    else:
        src = remapper.src_grid_info
        if 'type' not in src:
            raise ValueError(
                'None of the "src_from_*()" methods were called')
        src_descriptor = _get_descriptor(src)
        src_descriptor.format = remapper.format

    if remapper.dst_descriptor is not None:
        dst_descriptor = remapper.dst_descriptor
    else:
        dst = remapper.dst_grid_info

        if 'type' not in dst:
            raise ValueError(
                'None of the "dst_from_*()" methods were called')

        dst_descriptor = _get_descriptor(dst)
        dst_descriptor.format = remapper.format

    if remapper.map_filename is None:
        map_tool = remapper.map_tool
        prefixes = {'esmf': 'esmf', 'moab': 'mbtr'}
        suffixes = {
            'conserve': 'aave',
            'bilinear': 'bilin',
            'neareststod': 'neareststod',
        }
        suffix = f'{prefixes[map_tool]}{suffixes[remapper.method]}'

        remapper.map_filename = (
            f'map_{src_descriptor.mesh_name}_to_{dst_descriptor.mesh_name}'
            f'_{suffix}.nc'
        )
    map_tool = remapper.map_tool
    method = remapper.method
    if map_tool not in ['moab', 'esmf']:
        raise ValueError(
            f'Unexpected map_tool {map_tool}. Valid '
            f'values are "esmf" or "moab".'
        )

    if (isinstance(dst_descriptor, PointCollectionDescriptor) and
            method not in ['bilinear', 'neareststod']):
        raise ValueError(
            f'method {method} not supported for destination '
            f'grid of type PointCollectionDescriptor.'
        )

    if map_tool == 'moab' and method == 'neareststod':
        raise ValueError('method neareststod not supported by mbtempest.')

    remapper.src_descriptor = src_descriptor
    remapper.dst_descriptor = dst_descriptor
