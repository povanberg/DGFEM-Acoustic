#include <iostream>
#include <gmshUtils.h>
#include <Element.h>
#include <Face.h>
#include <algorithm>
#include <cassert>
#include "Mesh.h"
#include "gmsh.h"
#include "logger.h"

Mesh::Mesh(std::string name) :  name(name){

    this->dim = gmsh::model::getDimension();

    // Get all physical groups
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags, -1);
    // Map physical name to Tag/dim
    for(auto dimTag : dimTags) {
        gmsh::model::getPhysicalName(dimTag.first, dimTag.second, name);
        this->physDimTags[name] = dimTag;
    }

    // Get entities for domain
    gmsh::model::getEntitiesForPhysicalGroup(
            this->physDimTags["domain"].first,
            this->physDimTags["domain"].second,
            this->domainEntityTags);

    // Get entities for dirichelet
    gmsh::model::getEntitiesForPhysicalGroup(
            this->physDimTags["dirichelet"].first,
            this->physDimTags["dirichelet"].second,
            this->diricheletEntityTags);

    // Get elements in domain for all entities
    for(auto entityTag : domainEntityTags) {

        //----------------------------------------------------------
        // 1: Retrieve element tags and node tags
        //----------------------------------------------------------
        // Gmsh api call
        std::vector<int> elementTags;
        std::vector<int> nodeTags;
        int elementType = gmshUtils::getElementType(this->dim);
        gmsh::model::mesh::getElementsByType(elementType, elementTags, nodeTags, entityTag);
        // Object assignment
        std::vector<int>::const_iterator elIt = nodeTags.begin();
        int numNodes = (int) nodeTags.size() / elementTags.size();
        for(int i=0; i<elementTags.size(); ++i, elIt+=numNodes) {
            std::vector<int> elementNodes(elIt, elIt+numNodes);
            Element element(this->dim, elementTags[i], elementNodes);
            this->elements.push_back(element);
        }
        Log("Elements loaded");

        //----------------------------------------------------------
        // 2: Get Jacobians
        //----------------------------------------------------------
        // Gmsh api call
        std::vector<double> jacobians;
        std::vector<double> determinants;
        std::vector<double> points;
        int numGauss = 4; // hard coded for now, see gauss3
        gmsh::model::mesh::getJacobians(this->elements[0].getType(), "Gauss3",
                                        jacobians, determinants, points, entityTag);
        // Object assignment
        std::vector<double>::const_iterator jIt = jacobians.begin();
        std::vector<double>::const_iterator detIt = determinants.begin();
        std::vector<double>::const_iterator pIt = points.begin();
        for(int i=0; i<elementTags.size(); ++i, jIt+=9*numGauss, detIt+=numGauss, pIt+=3*numGauss) {
            std::vector<double> elementJacobian(jIt, jIt + 9*numGauss);
            std::vector<double> elementDet(detIt, detIt + numGauss);
            std::vector<double> elementPoints(pIt, pIt + 3*numGauss);
            this->elements[i].setJacobian(elementJacobian, elementDet, elementPoints, numGauss);
        }
        Log("Jacobian loaded");

        //----------------------------------------------------------
        // 3: Basis functions
        //----------------------------------------------------------
        // Gmsh api call
        int numComp;
        std::vector<double> basisFunctions;
        std::vector<double> gradBasisFunctions;
        std::vector<double> integrationPoints;
        gmsh::model::mesh::getBasisFunctions(this->elements[0].getType(), "Gauss3", "Lagrange",
                                             integrationPoints, numComp, basisFunctions);
        gmsh::model::mesh::getBasisFunctions(this->elements[0].getType(), "Gauss3", "GradLagrange",
                                             integrationPoints, numComp, gradBasisFunctions);
        // Object assignment
        for(int i=0; i<elementTags.size(); ++i)
            this->elements[i].setBasis(basisFunctions, gradBasisFunctions, integrationPoints);
        Log("Basis functions loaded");
    }

    //----------------------------------------------------------
    // 4: Retrieve face tags and nodes
    //----------------------------------------------------------
    // Gmsh api call
    std::string faceName = gmshUtils::getFaceFamilyName(this->dim-1, this->elements[0].getName());
    int faceNumNodes = gmshUtils::getFaceNumNodes(this->dim, this->elements[0].getOrder());
    int faceType = gmsh::model::mesh::getElementType(faceName, this->elements[0].getOrder());
    std::vector<int> faceTags;
    std::vector<int> faceNodeTags;
    if(this->dim < 3)
        gmsh::model::mesh::getElementEdgeNodes(elements[0].getType(), faceNodeTags, -1);
    else
        gmsh::model::mesh::getElementFaceNodes(elements[0].getType(), 3, faceNodeTags, -1); // Todo: 3 for trigs, 4 for quads.
    gmshUtils::makeUniqueInterfaces(faceNodeTags, faceNumNodes);
    int facesTag = gmsh::model::addDiscreteEntity(this->dim-1);
    gmsh::model::mesh::setElementsByType(this->dim-1, facesTag, faceType, {}, faceNodeTags);
    faceNodeTags.clear();
    gmsh::model::mesh::getElementsByType(faceType, faceTags, faceNodeTags, facesTag);
    // Object Assignement
    std::vector<int>::const_iterator fIt = faceNodeTags.begin();
    for(unsigned int i=0; i<faceTags.size(); ++i, fIt+=faceNumNodes) {
        std::vector<int> faceNodeTagsCurrent(fIt, fIt+faceNumNodes);
        Face face(faceTags[i], faceName, this->dim-1, faceNumNodes, faceType, faceNodeTagsCurrent);
        for(auto element = std::begin(this->elements); element!=std::end(this->elements); ++element) {
            bool hasFace = true;
            for (auto node : face.getNodeTags()) {
                if(!element->hasNode(node))
                    hasFace = false;
            }
            if(hasFace) {
                element->addFace(face);
                face.addElement(element->getTag());
            }
        }
    }
    Log("Faces loaded");

    //for (std::vector<int>::const_iterator i = this->faceTags.begin(); i != this->faceTags.end(); ++i)
    //    std::cout << *i << ' ';
    //std::cout << std::endl;
};


// Assemble Mesh mass matrix from element mass matrix
// For efficiency the mass matrix is sparse (block diagonal)
typedef Eigen::Triplet<double> T;
void Mesh::getMassMatrix(Eigen::SparseMatrix<double> M){
    int offset = 0;
    int MSize = this->elements.size()*this->elements[0].getNumNodes();
    std::vector<T> tripletList;
    Eigen::MatrixXd elMassMatrix;
    // Memory approximation reserve
    tripletList.reserve(MSize*MSize);
    // Filling
    M.resize(MSize, MSize);
    for(Element& el : this->elements){
        el.getMassMatrix(elMassMatrix);
        for(int i=0; i<elMassMatrix.rows(); ++i){
            for(int j=0; j<elMassMatrix.cols(); ++j){
                tripletList.push_back(T(offset+i,offset+j, elMassMatrix(i,j)));
            }
        }
        offset += el.getNumNodes();
    }
    // Triplets -> Sparse
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}


void Mesh::getStiffMatrix(Eigen::SparseMatrix<double> K, const Eigen::Vector3d a) {
    int offset = 0;
    int KSize = this->elements.size()*this->elements[0].getNumNodes();
    std::vector<T> tripletList;
    Eigen::MatrixXd elStiffMatrix;
    // Memory approximation reserve
    tripletList.reserve(KSize*KSize);
    // Filling
    K.resize(KSize, KSize);
    for(Element& el : this->elements){
        el.getStiffMatrix(elStiffMatrix, a);
        for(int i=0; i<elStiffMatrix.rows(); ++i){
            for(int j=0; j<elStiffMatrix.cols(); ++j){
                tripletList.push_back(T(offset+i,offset+j, elStiffMatrix(i,j)));
            }
        }
        offset += el.getNumNodes();
    }
    // Triplets -> Sparse
    K.setFromTriplets(tripletList.begin(), tripletList.end());
}