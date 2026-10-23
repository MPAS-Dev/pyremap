# PointCollectionDescriptor
```{index} single: Mesh Descriptors; PointCollectionDescriptor
```

A descriptor for a collection of points.

## Overview
The `PointCollectionDescriptor` class describes a set of points with latitude and longitude.

> **Note**: Point collections can only be used as destination meshes and are limited to the `bilinear` or `neareststod` remapping methods.

## Attributes
- `lat`: Latitude values of the points.
- `lon`: Longitude values of the points.

## Methods
- `to_scrip`: Converts the point collection to a SCRIP file.

## Example
```python
from pyremap import PointCollectionDescriptor

descriptor = PointCollectionDescriptor(lat, lon, "point_collection")
descriptor.to_scrip("points.scrip.nc")
```