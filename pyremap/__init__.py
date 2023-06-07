from pyremap.descriptor import (
    LatLon2DGridDescriptor,
    LatLonGridDescriptor,
    MpasEdgeMeshDescriptor,
    MpasMeshDescriptor,
    PointCollectionDescriptor,
    ProjectionGridDescriptor,
    get_lat_lon_descriptor,
)
from pyremap.polar import get_polar_descriptor, get_polar_descriptor_from_file
from pyremap.remapper import Remapper

__version_info__ = (1, 0, 0)
__version__ = '.'.join(str(vi) for vi in __version_info__)
