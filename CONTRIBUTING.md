# Contributing to GRO-DFTB

Thank you for your interest in contributing to GRO-DFTB.

## Reporting Issues

Please report bugs and feature requests via [GitHub Issues](https://github.com/ExaPsi/GRO-DFTB/issues).

When reporting a bug, include:
- GRO-DFTB version (`git describe --tags`)
- DFTB+ and GROMACS versions
- Operating system and compiler
- Minimal reproducing input files
- Full error output

## Contributing Code

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Write tests for your changes
4. Ensure all existing tests pass (`ctest --test-dir build`)
5. Submit a pull request

### Code Style

- C11 standard for `libgrodftb` (src/, include/)
- C++17 for GROMACS integration code
- Function names: `grodftb_` prefix for public API
- All public functions documented in headers

### Testing Requirements

- Every new function must have a corresponding test
- Finite-difference validation for any force-producing code path
- NVE energy conservation test for new embedding modes

## License

All contributions are licensed under LGPL-3.0-or-later, consistent with the project license.
