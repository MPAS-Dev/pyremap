import pyproj
import xarray as xr

from pyremap import (
    LatLon2DGridDescriptor,
    LatLonGridDescriptor,
    MpasCellMeshDescriptor,
    MpasEdgeMeshDescriptor,
    MpasVertexMeshDescriptor,
    PointCollectionDescriptor,
    ProjectionGridDescriptor,
    get_lat_lon_descriptor,
)


def _get_descriptor(info):
    """
    Get a mesh descriptor from the mesh info
    """
    grid_type = info['type']
    if grid_type == 'mpas':
        descriptor = _get_mpas_descriptor(info)
    elif grid_type == 'lon-lat':
        descriptor = _get_lon_lat_descriptor(info)
    elif grid_type == 'proj':
        descriptor = _get_proj_descriptor(info)
    elif grid_type == 'points':
        descriptor = _get_points_descriptor(info)
    else:
        raise ValueError(f'Unexpected grid type {grid_type}')

    return descriptor


def _get_mpas_descriptor(info):
    """
    Get an MPAS mesh descriptor from the given info
    """
    mesh_type = info['mpas_mesh_type']
    filename = info['filename']
    mesh_name = info['name']

    if mesh_type == 'cell':
        descriptor = MpasCellMeshDescriptor(
            filename=filename, mesh_name=mesh_name
        )
    elif mesh_type == 'vertex':
        descriptor = MpasVertexMeshDescriptor(
            filename=filename, mesh_name=mesh_name
        )
    elif mesh_type == 'edge':
        descriptor = MpasEdgeMeshDescriptor(
            filename=filename, mesh_name=mesh_name
        )
    else:
        raise ValueError(f'Unexpected MPAS mesh type {mesh_type}')

    return descriptor


def _get_lon_lat_descriptor(info):
    """
    Get a lon-lat descriptor from the given info
    """

    if 'dlat' in info and 'dlon' in info:
        lon_min = info['lon_min']
        lon_max = lon_min + 360.0
        descriptor = get_lat_lon_descriptor(
            dLon=info['dlon'],
            dLat=info['dlat'],
            lonMin=lon_min,
            lonMax=lon_max,
        )
    else:
        filename = info['filename']
        lon = info['lon']
        lat = info['lat']
        with xr.open_dataset(filename) as ds:
            lon_lat_1d = len(ds[lon].dims) == 1 and len(ds[lat].dims) == 1
            lon_lat_2d = len(ds[lon].dims) == 2 and len(ds[lat].dims) == 2
            if not lon_lat_1d and not lon_lat_2d:
                raise ValueError(
                    f'longitude and latitude coordinates {lon} '
                    f'and {lat} have unexpected sizes '
                    f'{len(ds[lon].dims)} and '
                    f'{len(ds[lat].dims)}.'
                )

        if lon_lat_1d:
            descriptor = LatLonGridDescriptor.read(
                filename=filename, lon_var_name=lon, lat_var_name=lat
            )
        else:
            descriptor = LatLon2DGridDescriptor.read(
                filename=filename, lon_var_name=lon, lat_var_name=lat
            )

    if 'name' in info:
        descriptor.mesh_name = info['name']

    return descriptor


def _get_proj_descriptor(info):
    """
    Get a ProjectionGridDescriptor from the given info
    """
    filename = info['filename']
    grid_name = info['name']
    x = info['x']
    y = info['y']
    if 'proj_attr' in info:
        with xr.open_dataset(filename) as ds:
            proj_str = ds.attrs[info['proj_attr']]
    else:
        proj_str = info['proj_str']

    proj = pyproj.Proj(proj_str)

    descriptor = ProjectionGridDescriptor.read(
        projection=proj,
        filename=filename,
        mesh_name=grid_name,
        x_var_name=x,
        y_var_name=y,
    )

    return descriptor


def _get_points_descriptor(info):
    """
    Get a PointCollectionDescriptor from the given info
    """
    filename = info['filename']
    collection_name = info['name']
    lon_var = info['lon']
    lat_var = info['lat']
    with xr.open_dataset(filename) as ds:
        lon = ds[lon_var].value
        lat = ds[lat_var].values
        unit_attr = lon.attrs['units'].lower()
        if 'deg' in unit_attr:
            units = 'degrees'
        elif 'rad' in unit_attr:
            units = 'radians'
        else:
            raise ValueError(f'Unexpected longitude unit unit {unit_attr}')

    descriptor = PointCollectionDescriptor(
        lons=lon, lats=lat, collectionName=collection_name, units=units
    )

    return descriptor
