#include <configParser.h>
#include <Mesh.h>

#ifndef DGALERKIN_SOLVER_H
#define DGALERKIN_SOLVER_H

namespace solver {
    /**
     * Solve using forward explicit scheme. O(h)
     *
     * @param u initial nodal solution vector
     * @param mesh
     * @param config
     */
    void forwardEuler(std::vector<std::vector<double>> &u, Mesh &mesh, Config config);

    /**
     * Solve using explicit Runge-Kutta integration method. O(h^4)
     *
     * @param u initial nodal solution vector
     * @param mesh
     * @param config
     */
    void rungeKutta(std::vector<std::vector<double>> &u, Mesh &mesh, Config config);
}

#endif