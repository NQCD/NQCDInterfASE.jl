# NQCDInterfaceASE

[![Build Status](https://github.com/alexsp32/NQCDInterfaceASE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alexsp32/NQCDInterfaceASE.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package contains various interfaces between NQCD packages and the Python-based [Atomic Simulation Environment](). 
Since `ase` is the standard choice of package for atomic structure modification, the goal is to provide an interface that allows for easy connection to `ase.calculators` for compatibility with various Python-based MLIP codes. 
In addition, structure conversion between the NQCD format and `ase` is provided for I/O to a larger number of file formats. 

This interface uses PythonCall to ensure a Python environment with `ase` installed is present. You can use your own Python version and environment by modifying CondaPkg / PythonCall settings as follows:

## Using your system Python and own environment

`JULIA_CONDAPKG_BACKEND`

## Compatibility with PyCall.jl

Various older Julia packages use the [PyCall.jl]() package to interface to Python. PyCall.jl and PythonCall.jl can coexist and be loaded at the same time, however they must be using
the same Python executable. To ensure this happens, make sure PyCall.jl is using the correct Python executable (and any package using it loads properly before NQCDInterfASE) is loaded. 

For PythonCall.jl and PythonCall to be loadable together, set `JULIA_PYTHONCALL_EXE=@PyCal` before launching Julia. 
