![pyremap logo](_static/logo.png)
(pyremap)=

# pyremap

Python remapping tools for climate and earth system models.

## Introduction
```{index} single: Introduction
```

`pyremap` is a Python library designed to facilitate the remapping
(interpolation) of data between different spatial representations, such as
grids and meshes. Remapping is a critical step in climate and Earth system
modeling, where data from various sources—such as models, observations, or
reanalysis products—often need to be compared, combined, or transformed to a
common spatial domain.

The library provides a high-level interface for defining source and
destination grids or meshes, generating mapping (weight) files, and applying
these mappings to remap data. It supports structured grids (e.g., regular
latitude-longitude grids) and unstructured meshes (e.g., MPAS meshes), making
it a versatile tool for a wide range of applications.

Key features of `pyremap` include:
- **Support for multiple grid types**: Structured grids, unstructured meshes,
  and projected grids.
- **Flexible remapping methods**: Bilinear, conservative, and nearest-neighbor
  interpolation.
- **Integration with external tools**: Leverages tools like ESMF (Earth System
  Modeling Framework) and MOAB for generating mapping files.
- **Ease of use**: High-level Python API for defining grids, building
  mappings, and remapping data.

Whether you are working with climate model output, observational data, or
custom spatial datasets, `pyremap` provides the tools you need to perform
accurate and efficient remapping.

## Documentation Overview
```{index} single: Documentation; Overview
```

This documentation is organized into the following sections:

```{toctree}
:maxdepth: 2

quick_start
mesh_descriptors/index
remapper/index
polar_projections/index
examples/index
developer_guide/index
```

- **Quick Start**: A brief guide to get you started with `pyremap`.
- **Mesh Descriptors**: Detailed descriptions of the different types of grids and meshes supported by `pyremap`.
- **Remapper**: An overview of the remapping process, including source and destination grid setup, mapping file generation, and data remapping.
- **Polar Projections**: Information on working with polar stereographic projections.
- **Examples**: Step-by-step examples demonstrating common use cases.
- **Developer Guide**: Resources for contributors, including API references, testing instructions, contribution guidelines, and coding standards.

# Indices and tables
```{index} single: Indices and tables
```

* {ref}`genindex`
