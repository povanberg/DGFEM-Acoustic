#include "Face.h"

Face::Face(int tag, std::string name, int dim, int numNodes, int type, std::vector<int> nodeTags){
    this->tag = tag;
    this->name = name;
    this->dim = dim;
    this->numNodes = numNodes;
    this->type = type;
    this->nodeTags = nodeTags;
}

Face &Face::addElement(int tag) {
    this->elementTags.push_back(tag);
    return *this;
}