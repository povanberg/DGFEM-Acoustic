//
// Created by pierre-olivier on 24/02/19.
//

#ifndef DGALERKIN_FACE_H
#define DGALERKIN_FACE_H


#include <string>
#include <vector>

class Face {
    private:
    public:
        Face(int tag, std::string name, int dim, int numNodes, int type, std::vector<int> nodeTags);
        int tag;
        std::string name;
        int numNodes;
        int type;
        int dim;
        std::vector<int> nodeTags;
        std::vector<int> elementTags;

        Face &addElement(int tag);
};


#endif //DGALERKIN_FACE_H
