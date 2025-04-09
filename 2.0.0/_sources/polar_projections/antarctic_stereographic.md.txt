# Antarctic Stereographic Projection

The Antarctic stereographic projection is widely used for mapping the Antarctic region. It is ideal for preserving geometric properties in polar areas.

## Usage in pyremap

The `get_antarctic_stereographic_projection()` method provides a common
`pyproj` projection that can be used in `pyremap` workflows and beyond.


```python
from pyremap.polar import get_antarctic_stereographic_projection

# Example: Create a polar descriptor for the Antarctic
antarctic_proj = get_antarctic_stereographic_projection()
antarctic_descriptor = get_polar_descriptor(
    lx=6000.0, ly=5000.0, dx=10.0, dy=10.0, projection=antarctic_proj)
antarctic_descriptor.to_scrip(f'{antarctic_descriptor.mesh_name}_scrip.nc')
```
