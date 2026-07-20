# Testing Instructions
```{index} single: Testing Instructions
```

This document provides detailed instructions for testing `pyremap`.

## Running Tests

1. Ensure all dependencies are installed. You can use `conda` to install any missing dependencies:
   ```bash
   conda install --file dev-spec.txt
   ```

2. Install pyremap in edit mode:
   ```bash
   python -m pip install --no-deps --no-build-isolation -e .
   ```

3. Run the test suite using `pytest`:
   ```bash
   pytest --pyargs tests
   ```

4. To run a specific test file, provide the file path:
   ```bash
   pytest path/to/test_file.py
   ```

5. For more verbose output, use the `-v` flag:
   ```bash
   pytest -v --pyargs tests
   ```

## Writing Tests

- Place all test files in the `pyremap/tests/` directory.
- Use descriptive names for test functions and files.
- Follow the `pytest` conventions for writing tests.

## Debugging Failures

- Use the `--pdb` flag to drop into the Python debugger on test failure:
  ```bash
  pytest --pdb --pyargs tests
  ```

- Review the stack trace provided by `pytest` to identify the source of the failure.

## Additional Resources

- Refer to the [pytest documentation](https://docs.pytest.org/en/latest/) for advanced usage and features.
