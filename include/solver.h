#ifndef DGALERKIN_SOLVER_H
#define DGALERKIN_SOLVER_H

#include <configParser.h>
#include <vector>
#include <Mesh.h>

namespace solver {

    void solveForwardEuler(Mesh &mesh, const Config config);
}

#endif //DGALERKIN_SOLVER_H
