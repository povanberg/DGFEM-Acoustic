# Discontinuous Galerkin Method for Acoustic Wave Propagation

This repository implements a discontinuous Galerkin finite element medthod (DGFEM) applied to the linearized Euler equations and the acoustic perturbation equations. The solver is based on [GMSH](http://gmsh.info/) library and supports a wide range of features:

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

## Authors

* Pierre-Olivier Vanberg
* Martin Lacroix
* Tom Servais
