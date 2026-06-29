# Developer Guide
```{index} single: Developer Guide
```

Welcome to the Developer Guide for `pyremap`. This section provides resources and guidelines for contributors.

## Coding Standards

- Follow [PEP 8](https://peps.python.org/pep-0008/) for Python code style.
- Use descriptive variable and function names.
- Write clear and concise docstrings for all public functions and classes.
- Ensure code is compatible with Python 3.8+.

## Testing

- Write unit tests for all new features and bug fixes.
- Use the `pytest` framework for testing.
- Place tests in the `pyremap/tests/` directory.
- Run tests locally before submitting changes:
  ```bash
  pytest --pyargs tests
  ```

## Contribution Guidelines

- Fork the repository and create a feature branch for your changes.
- Write clear and concise commit messages.
- Ensure all tests pass before submitting a pull request.
- Include documentation updates for any new features or changes.
- Submit a pull request with a detailed description of your changes.

## Additional Resources

```{toctree}
:maxdepth: 2

contribution_workflow
testing_instructions
generate_preview_docs
releasing
api
```
