# Contributing to VCFX

Thank you for your interest in contributing to VCFX! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

Please be respectful and considerate of others when contributing to this project. We aim to foster an inclusive and welcoming community.

## How to Contribute

### Reporting Bugs

If you find a bug in VCFX, please report it by creating an issue in our GitHub repository. When reporting a bug, please include:

- A clear, descriptive title
- A detailed description of the issue, including steps to reproduce
- The expected behavior
- The actual behavior observed
- Any relevant error messages or logs
- Your system information (OS, compiler version, etc.)

### Suggesting Enhancements

We welcome suggestions for new features or improvements to existing functionality. To suggest an enhancement:

1. Check if the enhancement has already been suggested or implemented
2. Create a new issue with a clear description of the enhancement
3. Explain why this enhancement would be useful to VCFX users

### Contributing Code

1. Fork the repository
2. Create a new branch for your feature or bug fix
3. Install the pre-commit hooks with `pre-commit install`
4. Write your code, following our coding standards. The hooks will automatically
   format changed C/C++ files.
5. Add tests for your changes
6. Ensure all tests pass
7. Update documentation as needed
8. Commit your changes with clear, descriptive commit messages
9. Submit a pull request

## Development Setup

### Code Formatting

We use `clang-format` to keep the C++ code style consistent. A basic configuration
is provided in `.clang-format` at the repository root. The `pre-commit` hooks
automatically apply `clang-format` to any changed C/C++ files. You can also run
the hook manually if needed:

```bash
pre-commit run clang-format --files path/to/changed_file.cpp
```
If running `pre-commit run --all-files` reports that the `clang-format`
hook modified files, you can apply the formatting to the entire
repository with:

```bash
pre-commit run clang-format --all-files
```

### Prerequisites

- CMake (version 3.10 or higher)
- C++11 compatible compiler (GCC, Clang, MSVC)
- Git

### Building for Development

```bash
git clone https://github.com/ieeta-pt/VCFX.git
cd VCFX
mkdir build && cd build
cmake ..
make
```

### Running Tests

After building the project, run the test suite from the `build` directory:

```bash
ctest --output-on-failure
```

You can still run all shell tests directly if needed:

```bash
bash ../tests/test_all.sh
```

## Coding Standards

- Use consistent indentation (4 spaces)
- Follow naming conventions:
  - Class names: CamelCase
  - Functions and methods: camelCase
  - Variables: snake_case
  - Constants: UPPER_SNAKE_CASE
- Write clear, descriptive comments
- Document public API methods
- Keep lines to a reasonable length (around 80-100 characters)
- Use descriptive variable and function names

## Documentation

- Update the appropriate documentation when changing functionality
- Document all public methods and classes
- Provide usage examples when relevant
- Use clear, concise language
- Check for spelling and grammar errors

## Pull Request Process

1. Update the README.md or relevant documentation with details of changes
2. Update the version number if applicable, following [Semantic Versioning](http://semver.org/)
3. The pull request will be merged once it has been reviewed and approved by a maintainer

## License

By contributing to VCFX, you agree that your contributions will be licensed under the project's MIT license. 
