#include <iostream>
#include "Face.h"
#include "gmsh.h"
#include <Eigen/Dense>

Face::Face(int tag, std::string name, int dim, int numNodes, int type, std::vector<int> &nodeTags){
    this->tag = tag;
    this->name = name;
    this->dim = dim;
    this->numNodes = numNodes;
    this->type = type;
    this->nodeTags = nodeTags;
    setNormal();
}

void Face::setNormal(){
    // Gmsh api call
    std::vector<double> coord1;
    std::vector<double> coord2;
    std::vector<double> paramCoords;
    gmsh::model::mesh::getNode(this->nodeTags.front(), coord1, paramCoords);
    gmsh::model::mesh::getNode(this->nodeTags.back(), coord2, paramCoords);
    // Compute normal
    switch(this->dim) {
        case 0:
            this->normal = {1, 0, 0};
            break;
        case 1:
            this->normal = {coord1[1]-coord2[1], coord2[0]-coord1[0], 0};
            this->normal.normalize();
            break;
        case 2:
            // Todo: -> gmsh
            break;
    }
}

// Add adjacent to the face
Face &Face::addElement(int tag) {
    this->elementTags.push_back(tag);
    return *this;
}

// Get Nodes tags
const std::vector<int> &Face::getNodeTags(){
    return this->nodeTags;
}

// Get Tag
const int &Face::getTag(){
    return this->tag;
}

// Get normal
const Eigen::Vector3d &Face::getNormal(){
    return this->normal;
}