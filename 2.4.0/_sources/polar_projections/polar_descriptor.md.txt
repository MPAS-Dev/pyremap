# Polar Descriptor

The polar descriptor provides metadata and configuration for polar (Arctic or Antarctic) projections in pyremap. It helps define the parameters and characteristics of these projections.

## Usage in pyremap

The polar descriptor can be used to define and manage projection parameters for polar regions.

```python
from pyremap.polar import get_polar_descriptor

# Example: Create a polar descriptor for the Arctic
arctic_descriptor = get_polar_descriptor(
    lx=6000.0, ly=5000.0, dx=10.0, dy=10.0, projection='arctic')
arctic_descriptor.to_scrip(f'{arctic_descriptor.mesh_name}_scrip.nc')
```
