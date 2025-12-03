# Pull Request Checklist

Thank you for contributing to VCFX! Before submitting your pull request, please confirm the following:

- [ ] I ran `pre-commit run --files <changed files>` to execute `ruff`, `flake8`, `mypy`, and `pytest`.
- [ ] All C++ tests pass via `ctest --output-on-failure` from the build directory.
- [ ] All Python tests pass via `pytest tests/python`.
- [ ] Documentation has been updated where applicable (e.g. `README.md`, `docs/*`).

Provide a brief description of the changes below.

