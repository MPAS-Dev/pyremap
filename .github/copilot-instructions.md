# pyremap

pyremap is a Python library for remapping climate and earth system model data between different grids and projections. It provides tools for creating mapping files and performing remapping operations using ESMF (Earth System Modeling Framework) and ncremap.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

### Bootstrap, build, and test the repository:
- Setup conda environment:
  ```bash
  conda config --add channels conda-forge
  conda config --set channel_priority strict
  conda create -y -n pyremap --file dev-spec.txt
  ```
  - NEVER CANCEL: Environment creation takes 15-20 minutes. Set timeout to 30+ minutes.
- Activate environment and install package:
  ```bash
  source /usr/share/miniconda/etc/profile.d/conda.sh  # if conda init wasn't run
  conda activate pyremap
  python -m pip install --no-deps --no-build-isolation -e .
  ```
  - Takes ~30 seconds to complete.
- Run tests:
  ```bash
  pytest --pyargs tests
  ```
  - NEVER CANCEL: Test suite takes ~15 seconds for 18 tests. Set timeout to 5+ minutes.
- Run tests with verbose output:
  ```bash
  pytest --pyargs tests -v
  ```

### Build documentation:
- Build docs:
  ```bash
  cd docs
  DOCS_VERSION=main make versioned-html
  ```
  - NEVER CANCEL: Documentation build takes ~4 minutes. Set timeout to 10+ minutes.
- Preview docs locally:
  ```bash
  cd _build/html
  python -m http.server 8000
  # Open http://0.0.0.0:8000/main/ in browser
  ```

### Linting and code quality:
- Run ruff linting check:
  ```bash
  ruff check .
  ```
  - Takes ~0.01 seconds.
- Auto-fix linting issues:
  ```bash
  ruff check --fix .
  ```
- Format code:
  ```bash
  ruff format .
  ```
  - Takes ~0.01 seconds.
- Run type checking:
  ```bash
  mypy --config=pyproject.toml pyremap
  ```

## Validation

- ALWAYS run through at least one complete end-to-end scenario after making changes.
- Test basic functionality by running:
  ```python
  from pyremap import Remapper, get_lat_lon_descriptor
  descriptor = get_lat_lon_descriptor(dlon=5.0, dlat=5.0)
  remapper = Remapper(ntasks=1, method='bilinear')
  remapper.src_descriptor = descriptor
  remapper.dst_descriptor = descriptor
  print("Basic functionality working")
  ```
- ALWAYS run `ruff check --fix .` and `ruff format .` before you are done or the CI (.github/workflows/build_workflow.yml) will fail.
- ALWAYS run the test suite with `pytest --pyargs tests` after making changes.

## Common tasks

The following are outputs from frequently run commands. Reference them instead of viewing, searching, or running bash commands to save time.

### Repository root structure
```
.github/         # GitHub workflows and configuration
ci/              # Conda build recipes
dev-spec.txt     # Development dependencies for conda
docs/            # Sphinx documentation
examples/        # Example scripts
pyremap/         # Main Python package
tests/           # Test suite
pyproject.toml   # Project configuration
README.md        # Project readme
```

### Key package modules
```
pyremap/
├── __init__.py           # Main package imports
├── descriptor/           # Grid descriptors (MPAS, lat-lon, projections)
├── remapper/             # Core remapping functionality
├── polar.py             # Polar projection utilities
├── utility.py           # Helper functions
└── version.py           # Version information
```

### Available external tools
- `ESMF_RegridWeightGen`: Available for creating mapping files with ESMF
- `ncremap`: Available for remapping NetCDF files
- `mbtempest`: NOT available (MOAB tool for mapping)

### Common pyremap workflows
1. **Create mapping file**: Use `Remapper.build_map()` to generate mapping weights
2. **Remap with ncremap**: Use `Remapper.ncremap()` to remap NetCDF files
3. **Remap with numpy**: Use `Remapper.remap_numpy()` for in-memory remapping
4. **Grid descriptors**: Use classes like `LatLonGridDescriptor`, `MpasCellMeshDescriptor`, etc.

### Dependencies and requirements
- Python >=3.9
- Core dependencies: dask, netcdf4, numpy, pyproj, scipy, xarray >=0.10.0
- External tools: ESMF (included), NCO (included for ncremap)
- Development tools: pytest, ruff, mypy, sphinx

### CI/CD information
- Tests run on Python 3.9-3.13
- Build timeout: 20 minutes
- Pre-commit hooks: ruff (linting/formatting), mypy (type checking)
- Documentation is auto-deployed from main branch

### Timing expectations
- **Environment setup**: 15-20 minutes (conda create)
- **Package installation**: ~30 seconds
- **Test suite**: ~15 seconds (18 tests)
- **Documentation build**: ~4 minutes
- **Linting/formatting**: <1 second
- **Basic import**: <1 second

### Manual validation scenarios
- Create a basic lat-lon grid descriptor: `get_lat_lon_descriptor(dlon=5.0, dlat=5.0)`
- Create a remapper instance: `Remapper(ntasks=1, method='bilinear')`
- Set source and destination grids
- Verify grid properties: `descriptor.mesh_name`, `descriptor.dims`, `descriptor.dim_sizes`

### Known limitations
- MOAB tools (mbtempest) are not available in the conda environment
- Network-dependent operations (intersphinx, pre-commit remote hooks) may fail in restricted environments
- Some examples require external MPAS grid files that are not included in the repository