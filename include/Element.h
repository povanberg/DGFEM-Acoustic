//
// Created by pierre-olivier on 24/02/19.
//

#ifndef DGALERKIN_ELEMENT_H
#define DGALERKIN_ELEMENT_H

#include <string>
#include <vector>
#include "Face.h"
#include "Eigen/Dense"

class Element {
    private:
        // Implements getter/setter, would be nice

    public:
        Element(int dim, int tag, std::vector<int> nodeTags);

        int dim;
        int tag;
        int type;
        int order;
        int numNodes;
        std::string name;
        std::vector<double> paramCoord;
        std::vector<int> nodeTags;
        std::vector<Face> faces;

        int numIntPoints;
        int numBasisFcts;
        std::vector<Eigen::Matrix3d> jacobian;
        std::vector<double> detJacobian;
        std::vector<Eigen::Vector3d> xPoints;    // x = physical coordinates
        std::vector<Eigen::Vector3d> uPoints;    // u = parametric coordinates
        std::vector<double> xbasisFct;
        std::vector<double> xgradBasisFct;  // see addBasis()
        std::vector<double> weights;

        bool hasNode(int tag);
        Element &addFace(Face face);
        Element &setJacobian(std::vector<double> &jacobian,
                             std::vector<double> &detJacobian,
                             std::vector<double> &xPoints,
                             int numIntPoints);
        Element &addBasis(std::vector<double> &ubasisFct,
                          std::vector<double> &ugradBasisFct,
                          std::vector<double> &uPoints);
};


#endif //DGALERKIN_ELEMENT_H
