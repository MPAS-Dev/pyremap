# Generating and Previewing Documentation
```{index} single: Documentation; Generating and Previewing
```

This guide explains how to build and preview the documentation for `pyremap`.

## Prerequisites

Ensure you have the following installed:
- [Miniforge](https://github.com/conda-forge/miniforge) or [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
- The required dependencies listed in `dev-spec.txt`:
  ```bash
  conda create -n pyremap_dev --file dev-spec.txt
  conda activate pyremap_dev
  python -m pip install --no-deps --no-build-isolation -e .
  ```

## Building the Documentation

To build the documentation, run the following command:
```bash
cd docs
DOCS_VERSION=main make versioned-html
```
This will generate the static site in the `_build/html/main` directory.

## Previewing the Documentation

To preview the documentation locally, open the `index.html` file in the `_build/html/main` directory with your browser or try:
```bash
cd _build/html
python -m http.server 8000
```
Then, open http://0.0.0.0:8000/main/ in your browser.

## Notes

- Make sure all documentation changes are committed before building or previewing.
- Verify the formatting and content in the generated HTML files before submitting changes.
