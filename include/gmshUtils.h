#include <gmsh.h>

#ifndef DGALERKIN_GMSHUTILS_H
#define DGALERKIN_GMSHUTILS_H

namespace gmshUtils {

    int getElementType(int dim);

    std::string getFaceFamilyName(const int dim, const std::string name);

    int getFaceNumNodes(int dim, int order);

    void makeUniqueInterfaces(std::vector<int> &nodes, const int n);

}

#endif
