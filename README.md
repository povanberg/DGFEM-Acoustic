# MATH0471-DG
This project is proposed as part of the course: multiphysics integrated computational project ([MATH 0471](http://www.montefiore.ulg.ac.be/~geuzaine/MATH0471/)) at Uliège university, spring 2019.

## Introduction

This  project  consists  in  studying  a  hyperbolic  system  of  equations  in  its  conservation form.   Spatial  discretization  is  performed  using  the  Discontinuous  Galerkin  (DG) method  and  Lagrange  nodal  basis  functions  on  unstructured  meshes. Time integration is performed using time marching methods. (Euler, Runge-Kutta)

## Compilation

### Libraries

First, make sure the following libraries are installed and the corresponding path are correctly exported.
```
Gmsh
Eigen
Lapack
Blas
OpenMP
```

### Build

```
git clone https://github.com/pvanberg/MATH0471-DG.git
cd MATH0471-DG
mkdir build && cd build
cmake .. && make -j4
```

### Run
Once the sources sucessfully build, you can start playing with the solver. It required two arguments: a mesh file created with Gmsh and a config file containing the sikver options.
```
cd bin
./dgalerkin mymesh.msh myconfig.conf
```

## Authors

* Pierre-Olivier Vanberg
* Martin Lacroix
* Tom Servais
