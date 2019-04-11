#include <configParser.h>
#include <Mesh.h>

#ifndef DGALERKIN_SOLVER_H
#define DGALERKIN_SOLVER_H

namespace solver {
    // Solve using forward explicit scheme. O(h)
    // u : initial solution vector    |   mesh   : ...
    // a : convection vector          |   config : ...
    void forwardEuler(std::vector<std::vector<double>> &u, std::vector<std::vector<std::vector<double>>> &a, Mesh &mesh, Config config);

    // Solve using explicit Runge-Kutta integration method. O(h^4)
    // u : initial solution vector    |   mesh   : ...
    // a : convection vector          |   config : ...
    void rungeKutta(std::vector<std::vector<double>> &u, Mesh &mesh, Config config);
}

#endif