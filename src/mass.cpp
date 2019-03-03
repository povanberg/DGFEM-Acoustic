#include "mass.h"
#include <vector>
#include <Mesh.h>
#include <iostream>
#include "Eigen/Dense"

namespace mass{

    void createM(Mesh &mesh){

        // Computation of the matrices N*N once, that are the same for each elements
        // One has a value of the Mass Matrix and of the set of shape functions of
        // the element at each of its Gauss points

        int nPts = mesh.elements[0].numIntPoints;
        int nFct = mesh.elements[0].numBasisFcts;
        std::vector<Eigen::MatrixXd> tempMlist;

        for(int k=0; k<nPts; ++k){
            Eigen::MatrixXd tempM(nFct,nFct);
            for(int i=0; i<nFct; ++i){
                for(int j=0; j<nFct; ++j){

                    tempM(i,j) = mesh.elements[0].xbasisFct[k](i)*mesh.elements[0].xbasisFct[k](j);
                }
            }
            tempMlist.push_back(tempM);
        }

        // Computation of M = sum of all w*N*N*detJ of the points in the element
        // Do this for each elements, w*N*N* is the same, but detJ is different

        for(int k=0; k<mesh.elements.size(); ++k){
             mesh.elements[k].M = mesh.elements[k].weights[0]*tempMlist[0]*mesh.elements[k].detJacobian[0];
            for(int i=1; i<nPts; ++i){
                mesh.elements[k].M += mesh.elements[k].weights[i]*tempMlist[i]*mesh.elements[k].detJacobian[i];
            }
        }
    }
}