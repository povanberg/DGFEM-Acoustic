#include "gmshUtils.h"
#include <gmsh.h>
#include <algorithm>

namespace gmshUtils {

    int getElementType(int dim) {
        std::vector<int> eleTypes;
        gmsh::model::mesh::getElementTypes(eleTypes, dim);
        if (eleTypes.size() != 1)
            gmsh::logger::write("Hybrid meshes not handled!", "error");
        return eleTypes[0];
    }

    std::string getFaceFamilyName(const int dim, const std::string name)
    {
        switch(dim)
        {
            case 0:
                return "point";
            case 1:
                return "line";
            case 2:
                if(name.substr(0,11) == "Tetrahedron")
                    return "triangle";
                else
                    return "quadrangle";
        }
    }

    int getFaceNumNodes(int dim, int order) {
        switch (dim)
        {
            case 1:
                return 1;
            case 2:
                return 2 + (order-1);
            case 3:
                return 3 + 3*(order-1);
        }
    }

    void makeUniqueInterfaces(std::vector<int> &nodes, const int n) {
        // Sort nodes per interfaces
        for(int i=0; i<nodes.size(); i+=n)
            std::sort(nodes.begin()+i, nodes.begin()+(i+n));

        // Remove identical interfaces
        for(std::vector<int>::iterator it = nodes.begin(); it != nodes.end();)
        {
            auto it_duplicate = std::search(it+n, nodes.end(), it, it+n);

            if(it_duplicate != nodes.end())
                nodes.erase(it_duplicate, it_duplicate+n);
            else
                it+=n;
        }
    }
}