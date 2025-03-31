### Smoothing with `expand_dist` and `expand_factor`

`pyremap` provides a unique capability for controlled smoothing when writing out the destination mesh descriptor. This is achieved through the use of two optional attributes: `expand_dist` and `expand_factor`. These attributes can be specified either as part of the destination descriptor or as attributes of the `remapper`, which will then pass them along to the destination descriptor.

#### Attributes Overview

- **`expand_dist`**:
  Specifies a distance in meters by which each grid cell is expanded. This can be provided as:
  - A single float value, which applies the same expansion distance to all grid cells.
  - A `numpy.ndarray` with the same size as the latitude and longitude coordinates of the mesh or grid, allowing for cell-specific expansion distances.

- **`expand_factor`**:
  Specifies a fraction of the cell size by which to expand each grid cell. This can also be provided as:
  - A single float value, applying a uniform expansion factor to all grid cells.
  - A `numpy.ndarray` with the same size as the latitude and longitude coordinates of the mesh or grid, enabling cell-specific expansion factors.

> **Note**: The default values are `expand_dist = 0.0` and `expand_factor = 1.0`, which result in no smoothing.

#### Usage and Behavior

These attributes are particularly useful for producing highly controlled smoothing as part of the remapping process. By expanding grid cells either by a fixed distance (`expand_dist`) or by a fraction of their size (`expand_factor`), users can fine-tune the smoothing effect to meet specific requirements.

#### Important Notes

1. **Combination with the `conserve` Method**:
   The `expand_dist` and `expand_factor` attributes must be used in combination with the `conserve` remapping method. This is because they leverage the area-weighted interpolation mechanism inherent to this method.

2. **Flexibility**:
   The ability to specify these attributes as either floats or `numpy.ndarray` fields provides flexibility for both uniform and non-uniform smoothing across the grid.

By using these attributes, `pyremap` enables advanced control over the remapping process, making it a powerful tool for applications requiring precise smoothing adjustments.
