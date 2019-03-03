//
// Created by pierre-olivier on 24/02/19.
//

#include <gmsh.h>
#include <gmshUtils.h>
#include <iostream>
#include "Element.h"
#include <algorithm>
#include <utils.h>
#include "Eigen/Dense"

Element::Element(int dim, int tag, std::vector<int> nodeTags) {

    // Global element informations (common to all elements, no hybrid mesh supported)
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

bool Element::hasNode(int tag){
    auto it = std::find(this->nodeTags.begin(), this->nodeTags.end(), tag);
    if(it != this->nodeTags.end())
        return true;
    else
        return false;
}

Element &Element::addFace(Face face){
    this->faces.push_back(face);
    return *this;
}

Element &Element::setJacobian(std::vector<double> &jacobian,
                              std::vector<double> &detJacobian,
                              std::vector<double> &xPoints,
                              int numIntPoints){

    for(int i=0; i<jacobian.size()/9; ++i){
        Eigen::Vector3d tempxPoints;
        Eigen::Matrix3d tempJacobian;
        for(int j=0; j<3; ++j){
            tempxPoints(j) = xPoints[i*3+j];
            for(int k=0; k<3; ++k){
                tempJacobian(j,k) = jacobian[i*9+j*3+k];
            }
        }
        this->xPoints.push_back(tempxPoints);
        this->jacobian.push_back(tempJacobian);
        this->invJacobian.push_back(tempJacobian.inverse());
    }
    this->detJacobian = detJacobian;
    this->numIntPoints = numIntPoints;
    this->numBasisFcts = numIntPoints;
    return *this;
}

Element &Element::addBasis(std::vector<double> &ubasisFct,
                           std::vector<double> &ugradBasisFct,
                           std::vector<double> &uPoints) {

    // Extract weigths
    for(int it = 0; it < uPoints.size(); it += 4)
        this->weights.push_back(uPoints[it+3]);
    // Remove weights from points
    for(int i=0; i<uPoints.size()/4; ++i){
        Eigen::Vector3d tempuPoints;
        Eigen::Vector3d tempFct;
        for(int j=0; j<3; ++j){
            tempuPoints(j) = uPoints[i*4+j];
            tempFct(j) = uPoints[i*3+j];
        }
        this->uPoints.push_back(tempuPoints);
        this->xbasisFct.push_back(tempFct);
    }

    // We have, df/du and want df/dx
    // solution: df/dx = df/du * du/dx
    // But jacobian is dx/du so inverse.
    int n = 3; // Can be optimized using the dimension, but as gmsh give 3D values...
    double dfdu, dfdx, dudx;
    for(int i=0; i<this->numIntPoints; ++i){

        // Inverse and store invJacobian (=du/dx)
        for(int j=0; j<this->numBasisFcts; ++j){
            for(int x=0; x<n; ++x) {
                dfdx = 0.;
                for(int u=0; u<n; ++u) {
                    dfdu = ugradBasisFct[i*n*numBasisFcts+j*numBasisFcts+u];
                    dudx = this->invJacobian[i](u,x);
                    dfdx += dfdu*dudx;
                }
                this->xgradBasisFct.push_back(dfdx);
            }
        }
    }

    return *this;
}