#include <gmsh.h>
#include <gmshUtils.h>
#include <iostream>
#include "Element.h"
#include <algorithm>
#include <logger.h>
#include <Eigen/Dense>

Element::Element(int dim, int tag, std::vector<int> &nodeTags) {

    this->dim = dim;
    this->tag = tag;
    this->nodeTags = nodeTags;
    this->type = gmshUtils::getElementType(this->dim);
    gmsh::model::mesh::getElementProperties(
            this->type,
            this->name,
            this->dim,
            this->order,
            this->numNodes,
            this->paramCoord);
}

// Add face to the element
Element &Element::addFace(Face face){
    this->faces.push_back(face);
    return *this;
}

// Set jacobian for the element
Element &Element::setJacobian(std::vector<double> &jacobian,
                              std::vector<double> &detJacobian,
                              std::vector<double> &xPoints,
                              int numIntPoints){

    this->numIntPoints = numIntPoints;
    this->numBasisFcts = this->numNodes;
    std::vector<double>::const_iterator gIt = jacobian.begin();
    for(unsigned int g=0; g<this->numIntPoints; ++g, gIt+=9){
        std::vector<double> gJacobian(gIt, gIt + 9);
        this->jacobian.push_back(Eigen::Map<Eigen::Matrix3d>(gJacobian.data()).transpose());
        this->invJacobian.push_back(this->jacobian[g].inverse());
    }
    this->detJacobian = Eigen::Map<Eigen::VectorXd>(detJacobian.data(), this->numIntPoints);
    this->xPoints = Eigen::Map<Eigen::MatrixXd>(xPoints.data(), 3, this->numIntPoints).transpose();

    return *this;
}

// Set basis functions for the element
Element &Element::setBasis(std::vector<double> &ubasisFct,
                           std::vector<double> &ugradBasisFct,
                           std::vector<double> &uPoints) {

    this->weights = Eigen::VectorXd(this->numIntPoints);
    this->uPoints = Eigen::MatrixXd(this->numIntPoints, 3);
    for(int g = 0; g < this->numIntPoints; ++g){
        this->uPoints(g, 0) = uPoints[g*4+0];
        this->uPoints(g, 1) = uPoints[g*4+1];
        this->uPoints(g, 2) = uPoints[g*4+2];
        this->weights(g) = uPoints[g*4+3];
    }

    Eigen::VectorXd ubasisFctEigen;
    Eigen::MatrixXd ugradBasisFctEigen;
    Eigen::MatrixXd xgradBasisFctEigen;
    for(int f = 0; f < this->numBasisFcts; ++f) {
        // Extract
        std::vector<double> ubFct(this->numIntPoints);
        std::vector<double> ubgradFct(this->numIntPoints*3);
        for(int g = 0; g < this->numIntPoints; ++g){
            ubFct[g] = ubasisFct[g*this->numBasisFcts+f];
            ubgradFct[g*3+0] = ugradBasisFct[g*this->numBasisFcts*3+f*3+0];
            ubgradFct[g*3+1] = ugradBasisFct[g*this->numBasisFcts*3+f*3+1];
            ubgradFct[g*3+2] = ugradBasisFct[g*this->numBasisFcts*3+f*3+2];
        }
        // To Eigen
        ubasisFctEigen = Eigen::Map<Eigen::VectorXd>(ubFct.data(), this->numIntPoints);
        ugradBasisFctEigen = Eigen::Map<Eigen::MatrixXd>(ubgradFct.data(), 3, this->numIntPoints).transpose();
        xgradBasisFctEigen = Eigen::MatrixXd(this->numIntPoints, 3);
        // df/du -> df/dx
        for(int g=0; g < this->numIntPoints; g++){
            xgradBasisFctEigen.row(g) = ugradBasisFctEigen.row(g)*this->invJacobian[g];
        }
        // Store
        this->basisFcts.push_back(ubasisFctEigen);
        this->gradBasisFcts.push_back(xgradBasisFctEigen);
    }
    return *this;
}

// Check if element contains the face
// with corresponding tag
bool Element::hasNode(const int tag){
    auto it = std::find(this->nodeTags.begin(), this->nodeTags.end(), tag);
    if(it != this->nodeTags.end())
        return true;
    else
        return false;
}

// Get element type
const int &Element::getTag(){
    return this->tag;
}

// Get element type
const int &Element::getType(){
    return this->type;
}

// Get element order
const int &Element::getOrder(){
    return this->order;
}

// Get element type name
const std::string &Element::getName(){
    return this->name;
}

// Get jacobian at gauss node g
const Eigen::Matrix3d &Element::getJacobian(const int g){
    return this->jacobian[g];
}

// Get inverse jacobian at gauss node g
const Eigen::Matrix3d &Element::getInvJacobian(const int g){
    return this->invJacobian[g];
};

// Get det jacobian for all g
const Eigen::VectorXd &Element::getDetJacobian(){
    return this->detJacobian;
}

// Get coordinates of gauss points
const Eigen::MatrixXd &Element::getXPoints(){
    return this->xPoints;
};

// Get parametric coordinates of gauss points
const Eigen::MatrixXd &Element::getUPoints(){
    return this->uPoints;
}

// Get weigths for integration points
const Eigen::VectorXd &Element::getWeigths(){
    return this->weights;
}

// Get Basis functions 'f' in numBasisFcts
const Eigen::VectorXd &Element::getBasisFcts(const int f){
    return this->basisFcts[f];
}

// Get Gradient Basis functions 'f' in numBasisFcts
const Eigen::MatrixXd &Element::getGradBasisFcts(const int f){
    return this->gradBasisFcts[f];
}

// Get number of nodes
const int &Element::getNumNodes(){
    return this->numNodes;
}

// Get element mass matrix
void Element::getMassMatrix(Eigen::MatrixXd &massMatrix){
    massMatrix.resize(this->numBasisFcts, this->numBasisFcts);
    for(unsigned int i=0; i<this->numBasisFcts; ++i){
        for(unsigned int j=0; j<this->numBasisFcts; ++j){
            massMatrix(i,j) = weights.cwiseProduct(basisFcts[i])
                                     .cwiseProduct(basisFcts[j])
                                     .dot(detJacobian);
        }
    }
}

// Get element mass matrix
void Element::getStiffMatrix(Eigen::MatrixXd &stiffMatrix, const Eigen::Vector3d a){
    stiffMatrix.resize(this->numBasisFcts, this->numBasisFcts);
    for(unsigned int i=0; i<this->numBasisFcts; ++i){
        for(unsigned int j=0; j<this->numBasisFcts; ++j){
            stiffMatrix(i,j) = weights.cwiseProduct(gradBasisFcts[i]*a)
                                      .cwiseProduct(basisFcts[j])
                                      .dot(detJacobian);
        }
    }
}