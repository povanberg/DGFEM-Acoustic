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
    setBasisFcts();
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

void Face::setBasisFcts() {
    // Gmsh api call
    int numComp;
    std::vector<double> basisFcts;
    std::vector<double> uPoints;
    gmsh::model::mesh::getBasisFunctions(this->type, "Gauss3", "Lagrange",
                                         uPoints, numComp, basisFcts);
    // Eigen conversion
    this->numBasisFcts = this->numNodes;
    this->numIntPoints = (int) basisFcts.size() / this->numBasisFcts;
    this->basisFcts = Eigen::Map<Eigen::MatrixXd>(basisFcts.data(), this->numBasisFcts, this->numIntPoints);

    this->weights = Eigen::VectorXd(this->numIntPoints);
    this->uPoints = Eigen::MatrixXd(this->numIntPoints, 3);
    for(int g = 0; g < this->numIntPoints; ++g){
        this->uPoints(g, 0) = uPoints[g*4+0];
        this->uPoints(g, 1) = uPoints[g*4+1];
        this->uPoints(g, 2) = uPoints[g*4+2];
        this->weights(g) = uPoints[g*4+3];
    }
}

void Face::setJacobian(std::vector<double> &jacobian,
                       std::vector<double> &detJacobian,
                       std::vector<double> &xPoints){

    this->numIntPoints = numIntPoints;
    this->numBasisFcts = this->numNodes;
    std::vector<double>::const_iterator gIt = jacobian.begin();
    for(unsigned int g=0; g<this->numIntPoints; ++g, gIt+=9){
        std::vector<double> gJacobian(gIt, gIt + 9);
        this->jacobian.push_back(Eigen::Map<Eigen::Matrix3d>(gJacobian.data()).transpose());
    }
    this->detJacobian = Eigen::Map<Eigen::VectorXd>(detJacobian.data(), this->numIntPoints);
    this->xPoints = Eigen::Map<Eigen::MatrixXd>(xPoints.data(), 3, this->numIntPoints).transpose();
}

// Add adjacent to the face
Face &Face::addElement(int tag) {
    this->elementTags.push_back(tag);
    if(this->elementTags.size() > 1)
        this->boundary = false;
    else
        this->boundary = true;
    return *this;
}

bool Face::hasNode(const int tag){
    auto it = std::find(this->nodeTags.begin(), this->nodeTags.end(), tag);
    if(it != this->nodeTags.end())
        return true;
    else
        return false;
}

int Face::getSecondElement(const int tag1){
    if(tag1 == this->elementTags[0])
        return this->elementTags[1];
    else
        return this->elementTags[0];
}

double Face::getFluxInt(const Eigen::Vector3d &a){
    return a.dot(this->normal)*weights.cwiseProduct(basisFcts.row(0).transpose())
                                      .cwiseProduct(basisFcts.row(1).transpose())
                                      .dot(detJacobian);
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