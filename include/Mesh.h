#ifndef DGALERKIN_MESH_H
#define DGALERKIN_MESH_H

#include <gmsh.h>
#include <map>
#include "Element.h"

class Mesh {

    private:
        std::string name;

    public:
        // Constructor
        Mesh(std:: string name);
        // Internal variables
        // Todo: public variables anti-pattern -> getter and setter
        int dim;                                                    // Mesh dimension

        std::map<std::string, std::pair<int,int>> physDimTags;      // Associate physical name to tag/dim
        std::vector<int> domainEntityTags;                          // List of all geometrical entities in domain
        std::vector<int> diricheletEntityTags;                      // List of all geometrical entities in dirichelet

        std::vector<Element> elements;

        int faceDim;
        std::string faceName;
        int faceNumNodes;
        int faceType;
        int facesTag;
        std::vector<int> faceTags;
        std::vector<int> faceNodeTags;

};


#endif //DGALERKIN_MESH_H
