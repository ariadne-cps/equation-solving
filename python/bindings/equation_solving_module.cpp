#include "pybind11.hpp"
void equation_solving_submodule(pybind11::module&);
void equation_solving_module(pybind11::module& module) { equation_solving_submodule(module); }
