# DGFEM for Acoustic Wave Propagation 

[![Build Status](https://travis-ci.org/pvanberg/MATH0471-DG.svg?branch=master)](https://travis-ci.org/pvanberg/MATH0471-DG)   [![Maintenance](https://img.shields.io/badge/Version-1.1.0-e67e22.svg)](https://github.com/pvanberg/MATH0471-DG/releases/tag/v1.1.0) [![Maintenance](https://img.shields.io/badge/c++-14%20|%2017%20|%2020-27ae60.svg)](https://github.com/pvanberg/MATH0471-DG/releases/tag/v1.0.0) 

This repository implements a discontinuous Galerkin finite element method (DGFEM) applied to the linearized Euler equations and the acoustic perturbation equations. The solver is based on [GMSH](http://gmsh.info/) library and supports a wide range of features:

- 1D, 2D, 3D problems
- 4-th order Runge-Kutta
- High order elements
- Absorbing and reflecting boundaries
- Complex geometry and unstructured grid

For more information, a detailled report is available here(soon).

## Getting Started

### Prerequisites

First, make sure the following libraries are installed. If you are running a linux distribution (ubuntu, debian, ...), an installation [script](https://github.com/pvanberg/MATH0471-DG/blob/master/build.sh) is provided. 

```
Gmsh
Eigen
Lapack
Blas
OpenMP
```

### Installing

```
git clone https://github.com/pvanberg/MATH0471-DG.git
cd MATH0471-DG
mkdir build && cd build
cmake .. && make -j4
```

## Running the tests
Once the sources sucessfully build, you can start using with the solver. It required two arguments: a mesh file created with Gmsh and a config file containing the solver options. Examples of mesh files and config files are given [here](https://github.com/pvanberg/MATH0471-DG/tree/master/doc).

```
cd bin
./dgalerkin mymesh.msh myconfig.conf
```

### Minimal working example

2D propagation of an Gaussian initial condition over a square.

```
./dgalerkin ../../doc/2d/square.msh ../../doc/config/config.conf 
```

## Authors

* Pierre-Olivier Vanberg
* Martin Lacroix
* Tom Servais
