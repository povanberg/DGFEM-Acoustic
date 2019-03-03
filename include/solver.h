#include <configParser.h>
#include <vector>
#include <Mesh.h>

#ifndef DGALERKIN_SOLVER_H
#define DGALERKIN_SOLVER_H

namespace solver {

    void solveForwardEuler(Mesh &mesh, const Config config);
}

#endif
