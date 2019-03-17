#include <configParser.h>
#include <Mesh.h>

#ifndef DGALERKIN_SOLVER_H
#define DGALERKIN_SOLVER_H

namespace solver {

    void solveForwardEuler(Mesh &mesh, Config config);
}

#endif