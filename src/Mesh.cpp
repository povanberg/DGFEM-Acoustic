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
        std::vector<int> elementTags;
        std::vector<int> nodeTags;
        int elementType = gmshUtils::getElementType(this->dim);
        gmsh::model::mesh::getElementsByType(elementType, elementTags, nodeTags, entityTag);
        // Split vector O(n), DG treats elements independently
        // Should not impact the performances.
        // O(nlog(n)) achievable if ordering.
        std::vector<int>::const_iterator first;
        std::vector<int>::const_iterator last;
        int numNodes = (int) nodeTags.size() / elementTags.size();
        for(int i=0; i<elementTags.size(); ++i) {
            first = nodeTags.begin() + i*numNodes;
            last = nodeTags.begin() + (i+1)*numNodes;
            std::vector<int> elementNodes(first, last);
            Element element(this->dim, elementTags[i], elementNodes);
            this->elements.push_back(element);
        }

        // Jacobians
        std::vector<double> jacobians;
        std::vector<double> determinants;
        std::vector<double> points;
        int numGauss = 3; // hard coded for now, see gauss3
        gmsh::model::mesh::getJacobians(
                this->elements[0].type,      // All element same type
                "Gauss2",                    // For now assume Gauss3
                jacobians, determinants, points, entityTag);
        std::vector<double>::const_iterator firstD;
        std::vector<double>::const_iterator lastD;
        for(int i=0; i<elementTags.size(); ++i) {
            firstD = jacobians.begin() + i*(9*numGauss);
            lastD = jacobians.begin() + (i+1)*(9*numGauss);
            std::vector<double> elementJacobian(firstD, lastD);
            firstD = determinants.begin() + i*numGauss;
            lastD = determinants.begin() + (i+1)*numGauss;
            std::vector<double> elementDet(firstD, lastD);
            firstD = points.begin() + i*numGauss;
            lastD = points.begin() + (i+1)*numGauss;
            std::vector<double> elementPoints(firstD, lastD);
            this->elements[i].setJacobian(elementJacobian, elementDet, elementPoints, numGauss);
        }

        // Basis function
        std::vector<double> basisFunctions;
        std::vector<double> gradBasisFunctions;
        std::vector<double> integrationPoints;
        int numBasisFunction;
        gmsh::model::mesh::getBasisFunctions(this->elements[0].type, "Gauss2", "Lagrange",
                                             integrationPoints, numBasisFunction, basisFunctions);
        gmsh::model::mesh::getBasisFunctions(this->elements[0].type, "Gauss2", "GradLagrange",
                                             integrationPoints, numBasisFunction, gradBasisFunctions);
        assert(numBasisFunction==this->elements[0].numNodes);
        for(int i=0; i<elementTags.size(); ++i) {
            this->elements[i].addBasis(basisFunctions, gradBasisFunctions, integrationPoints);
        }
    }

    // Create faces
    std::string faceName = gmshUtils::getFaceFamilyName(this->dim-1, this->elements[0].name);
    int faceNumNodes = gmshUtils::getFaceNumNodes(this->dim, this->elements[0].order);
    int faceType = gmsh::model::mesh::getElementType(faceName, this->elements[0].order);
    std::vector<int> faceTags;
    std::vector<int> faceNodeTags;
    if(this->dim < 3)
        gmsh::model::mesh::getElementEdgeNodes(elements[0].type, faceNodeTags, -1);
    else
        gmsh::model::mesh::getElementFaceNodes(elements[0].type, 3, faceNodeTags, -1); // Todo: 3 for trigs, 4 for quads.
    // Gmsh return duplicated faces for internal elements -> make them unique.
    gmshUtils::makeUniqueInterfaces(faceNodeTags, faceNumNodes);
    // Create new entity
    int facesTag = gmsh::model::addDiscreteEntity(this->dim-1);
    gmsh::model::mesh::setElementsByType(this->dim-1, facesTag, faceType, {}, faceNodeTags);
    // Let Gmsh autogenerate the Tags, then retrieve all faces
    faceNodeTags.clear();
    gmsh::model::mesh::getElementsByType(faceType, faceTags, faceNodeTags, facesTag);
    // Add faces to corresponding Elements
    std::vector<int>::const_iterator first;
    std::vector<int>::const_iterator last;
    for(unsigned int i=0; i<faceTags.size(); ++i) {
        first = faceNodeTags.begin() + i*faceNumNodes;
        last = faceNodeTags.begin() + (i+1)*faceNumNodes;
        std::vector<int> faceNodeTagsCurrent(first, last);
        // Create the face
        Face face(faceTags[i], faceName, this->dim-1, faceNumNodes, faceType, faceNodeTagsCurrent);
        // Get Elements which contains the face
        for(auto element = std::begin(this->elements); element!=std::end(this->elements); ++element) {
            bool hasFace = true;
            for (auto node : face.nodeTags) {
                if(!element->hasNode(node))
                    hasFace = false;
            }
            if(hasFace) {
                // Each element has face instances
                element->addFace(face);
                // Each surface keep track of its neighbors
                face.addElement(element->tag);
            }
        }
    }

    //for (std::vector<int>::const_iterator i = this->faceTags.begin(); i != this->faceTags.end(); ++i)
    //    std::cout << *i << ' ';
    //std::cout << std::endl;
};
