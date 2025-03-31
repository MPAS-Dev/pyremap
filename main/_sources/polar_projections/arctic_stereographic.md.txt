# Arctic Stereographic Projection

The Arctic stereographic projection is commonly used for mapping the Arctic region. It preserves angles and shapes over small areas, making it suitable for polar studies.

## Usage in pyremap

The `get_arctic_stereographic_projection()` method provides a common `pyproj`
projection that can be used in `pyremap` workflows and beyond.

```python
from pyremap.polar import get_arctic_stereographic_projection

# Example: Create a polar descriptor for the Arctic
arctic_proj = get_arctic_stereographic_projection()
arctic_descriptor = get_polar_descriptor(
    lx=6000.0, ly=5000.0, dx=10.0, dy=10.0, projection=arctic_proj)
arctic_descriptor.to_scrip(f'{arctic_descriptor.mesh_name}_scrip.nc')

```
