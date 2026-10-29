# Ariadne Equation Solving

[![macOS Status](https://github.com/ariadne-cps/equation-solving/actions/workflows/macos.yml/badge.svg)](https://github.com/ariadne-cps/equation-solving/actions/workflows/macos.yml)
[![Ubuntu Status](https://github.com/ariadne-cps/equation-solving/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ariadne-cps/equation-solving/actions/workflows/ubuntu.yml)
[![Windows Status](https://github.com/ariadne-cps/equation-solving/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/equation-solving/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/equation-solving/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/equation-solving/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/equation-solving/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/equation-solving/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/equation-solving/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/equation-solving)

Ariadne Equation Solving is a standalone C++20 library for validated solution of systems of equations, fixed-point problems, and implicit functions. It provides Interval Newton, Krawczyk and Factored Krawczyk methods.

## Layout

Public headers live under `include/equation-solving`, implementations directly in `src`, C++ tests in `tests`, Python bindings and tests in `python`, and standalone C++ and Python tutorials in `tutorials` and `python/tutorials`. The umbrella header is `include/ariadne-equation-solving.hpp`.

## Dependencies

- [ariadne-cps/function](https://github.com/ariadne-cps/function), included as a Git submodule.
- [ariadne-cps/configuration](https://github.com/ariadne-cps/configuration), used for shared build configuration.

## Build

```sh
git clone --recurse-submodules https://github.com/ariadne-cps/equation-solving.git
cd equation-solving
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel --target tests
ctest --test-dir build --output-on-failure
```

## Python bindings

A compatible Python 3.11+ development environment enables the standalone `pyariadne` extension with bindings for the equation-solving algorithms.

## Tutorials

- C++: `tutorials/tutorial_equation_solving/tutorial_equation_solving.cpp`.
- Python: `python/tutorials/tutorial_equation_solving.py`.

## License

GNU General Public License v3.0.
