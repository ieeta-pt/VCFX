# Publishing VCFX Python Package to PyPI

This guide explains how to publish the VCFX Python package to PyPI (Python Package Index).

## Overview

The VCFX Python package (`vcfx`) provides Python bindings and convenient wrappers for the VCFX toolkit. The package can be published to:

1. **TestPyPI** - For testing the publishing process
2. **PyPI** - For official releases

## Prerequisites

Before publishing, ensure:

1. You have the necessary permissions on GitHub to create releases or run workflows
2. The version number has been updated appropriately in the codebase
3. All tests are passing in the CI/CD pipeline

## Publishing Methods

### Method 1: Automatic Publishing via GitHub Release (Recommended)

This is the recommended method for official releases:

1. **Create a new release on GitHub:**
   - Go to the repository's [Releases page](https://github.com/ieeta-pt/VCFX/releases)
   - Click "Draft a new release"
   - Create a new tag with the format `v1.0.3` (following semantic versioning)
   - Write release notes
   - Publish the release

2. **The automated workflow will:**
   - Build the source distribution and wheels
   - Run package validation checks
   - Publish to PyPI automatically
   - Also trigger Docker image publishing

### Method 2: Manual Publishing via GitHub Actions

For testing or manual control:

1. Go to the [Actions tab](https://github.com/ieeta-pt/VCFX/actions)
2. Select "Publish Python Package" workflow
3. Click "Run workflow"
4. Choose either `testpypi` or `pypi` as the target
5. Click "Run workflow" button

### Method 3: Local Publishing

For development or emergency releases:

```bash
# Install dependencies
pip install build twine

# Clean previous builds
cd python
rm -rf dist/ build/ *.egg-info

# Build the package
python -m build

# Check the package
twine check dist/*

# Upload to TestPyPI (for testing)
twine upload --repository testpypi dist/*

# Upload to PyPI (for production)
twine upload dist/*
```

## Version Management

The package version is automatically extracted from the CMakeLists.txt file. Before publishing:

1. Update the version in `CMakeLists.txt`:
   ```cmake
   project(VCFX VERSION 1.0.3)
   ```

2. Commit the version change:
   ```bash
   git add CMakeLists.txt
   git commit -m "Bump version to 1.0.3"
   git push
   ```

## Testing the Published Package

After publishing to TestPyPI:

```bash
# Install from TestPyPI
pip install -i https://test.pypi.org/simple/ vcfx

# Test the package
python -c "import vcfx; print(vcfx.get_version())"
```

After publishing to PyPI:

```bash
# Install from PyPI
pip install vcfx

# Test the package
python -c "import vcfx; print(vcfx.get_version())"
```

## Troubleshooting

### Build Failures

If the package fails to build:
- Check that CMakeLists.txt exists in the root directory
- Ensure all Python files have proper formatting (run pre-commit)
- Verify that the version extraction script works

### Publishing Failures

If publishing fails:
- Check that GitHub Actions has the necessary permissions
- Verify that the package name isn't already taken on PyPI
- Ensure the version number hasn't been used before

### Missing C++ Extension

The package can work without the C++ extension (pure Python mode). If CMake isn't available during installation, only the Python wrappers will be installed.

## Security Notes

- The GitHub Actions workflow uses PyPI's trusted publishing (OIDC)
- No API tokens are stored in the repository
- Ensure you have the appropriate GitHub environments configured:
  - `testpypi` environment for TestPyPI publishing
  - `pypi` environment for PyPI publishing

## Package Contents

The published package includes:

- Python bindings for VCFX C++ library (when built with CMake)
- Tool wrappers for all VCFX command-line tools
- Type hints and py.typed marker
- Comprehensive convenience functions

## Next Steps

After publishing:

1. Update the documentation to reflect the new version
2. Create a GitHub release if not already done
3. Announce the release to users
4. Update any dependent projects 

For additional guidance or questions, consult the maintainers via GitHub issues.
