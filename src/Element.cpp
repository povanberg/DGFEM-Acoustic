#include <gmsh.h>
#include <gmshUtils.h>
#include <iostream>
#include "Element.h"
#include <algorithm>
#include <logger.h>
#include <Eigen/Dense>

Element::Element(int dim, int tag, std::vector<int> &nodeTags, Eigen::Vector3d &barycenter) {

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
    this->u.resize(this->numNodes);
    this->barycenter = barycenter;
}

// Add face to the element
Element &Element::addFace(Face face){
    this->faces.push_back(face);
    return *this;
}

void Element::setNormals(){

    for(Face& fa : this->faces){
        Eigen::Vector3d normal;
        Eigen::Vector3d nTest;
        std::vector<double> coord1;
        std::vector<double> coord2;
        std::vector<double> paramCoords;
        gmsh::model::mesh::getNode(fa.getNodeTags().front(), coord1, paramCoords);
        gmsh::model::mesh::getNode(fa.getNodeTags().back(), coord2, paramCoords);

        switch(fa.getDim()){
            case 0:
                // Not done yet
                break;
            case 1:
                normal = {coord1[1]-coord2[1], coord2[0]-coord1[0], 0};
                nTest = {coord1[0]-barycenter[0], coord1[1]-barycenter[1], 0};
                if(normal.dot(nTest)<0){
                    normal = -normal;
                }
                normal.normalize();
                break;
            case 2:
                // Not done yet
                break;
        }
        this->normals.push_back(normal);
    }
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
        this->invJacobian.push_back(this->jacobian[g].inverse().transpose()); // TRANSPOSE INV JACOBIAN
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

// Get list of node tags
std::vector<int> &Element::getNodeTags(){
    return this->nodeTags;
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
void Element::getStiffMatrix(Eigen::MatrixXd &stiffMatrix, const Eigen::Vector3d &a){
    stiffMatrix.resize(this->numBasisFcts, this->numBasisFcts);
    for(unsigned int i=0; i<this->numBasisFcts; ++i){
        for(unsigned int j=0; j<this->numBasisFcts; ++j){
            stiffMatrix(i,j) = weights.cwiseProduct(gradBasisFcts[i]*a)
                                      .cwiseProduct(basisFcts[j])
                                      .dot(detJacobian);
        }
    }
}
/*
void Element::getFlux(Eigen::VectorXd &Flux, const Eigen::Vector3d &a, std::vector<Element> &elements){
    // Independent of numerical flux
    Eigen::VectorXd numDataIn;
    Eigen::VectorXd numDataOut;
    Flux = Eigen::VectorXd::Zero(this->numBasisFcts);
    for(unsigned int i=0; i<this->numBasisFcts; ++i) {
        for (unsigned int j = 0; j < this->numBasisFcts; ++j) {
            for(Face &face: this->faces){
                if(face.boundary){
                    // Not flux at boundaries
                    Flux(i) += 0.0;
                }
                else {
                    // If face is not a boundary
                    // We get is element other side of face.
                    int elOutTag = face.getSecondElement(this->tag);
                    int idElOut;
                    for(int i=0; i< elements.size(); ++i){
                        if(elements[i].getTag() == tag)
                            idElOut=  i;
                    }
                    Element elOut = elements[idElOut];

                    // If Node 'i' is in Surface
                    if(face.hasNode(this->nodeTags[i])){
                        // Get indice of out solution
                        for(int k=0; k<elOut.getNumNodes(); ++k) {
                            if (this->getNodeTags()[i] == elOut.getNodeTags()[k]) {
                                this->getData(numDataIn);
                                elOut.getData(numDataOut);
                                Flux(i) += ((numDataIn(i) + numDataOut(k)) / 2.)*face.getFluxInt(a);
                            }
                        }
                    }
                    else{
                        // Shape function is zero on surface
                        Flux(i) += 0.0;
                    }
                }

            }
        }
    }
}
*/
void Element::getFlux(Eigen::VectorXd &Flux, const Eigen::Vector3d &a, std::vector<Element> &elements, Eigen::VectorXd &u){
    Flux = Eigen::VectorXd::Zero(this->numNodes);
    for(int i=0; i<this->faces.size(); ++i){
        Eigen::Vector2d faceFlux;
        Eigen::Vector2d uFace;
        uFace.setZero();
        if(i=0){
            uFace(0) += u(0);
            uFace(1) += u(1);
            this->faces[i].getFluxInt(faceFlux, uFace, a);
            Flux(0) += faceFlux(0);
            Flux(1) += faceFlux(1);
        }
        if(i=1){
            uFace(0) += u(1);
            uFace(1) += u(2);
            this->faces[i].getFluxInt(faceFlux, uFace, a);
            Flux(1) += faceFlux(0);
            Flux(2) += faceFlux(1);
        }
        if(i=2){
            uFace(0) = u(0);
            uFace(1) = u(2);
            this->faces[i].getFluxInt(faceFlux, uFace, a);
            Flux(0) += faceFlux(0);
            Flux(2) += faceFlux(1);
        }
    }
}

void Element::getData(Eigen::VectorXd &data){
    data = this->u;
}

void Element::setData(Eigen::VectorXd &data){
    this->u = data;
}