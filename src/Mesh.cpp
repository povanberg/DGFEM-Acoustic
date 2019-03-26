#include <string>
#include <gmsh.h>
#include <assert.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <omp.h>

#include "configParser.h"
#include "Mesh.h"
#include "utils.h"

// Mesh constructor: load the mesh parameters
// from Gmsh api and create the elements mapping.
Mesh::Mesh(std::string name, Config config) :  name(name), config(config) {

    //---------------------------------------------------------------------
    // Elements
    //---------------------------------------------------------------------
    m_elDim = gmsh::model::getDimension();
    gmsh::model::mesh::getElementTypes(m_elType, m_elDim);
    gmsh::model::mesh::getElementProperties(m_elType[0], m_elName, m_elDim,
                                            m_elOrder, m_elNumNodes, m_elParamCoord);
    gmsh::model::mesh::getElementsByType(m_elType[0], m_elTags, m_elNodeTags);

    // Gauss quadrature performs exact integration of
    // polynomials of order n with p=(n+1)/2 integration points.
    m_elNum = m_elTags.size();
    m_elIntType = "Gauss" + std::to_string(2*m_elOrder);
    gmsh::model::mesh::getJacobians(m_elType[0], m_elIntType, m_elJacobians,
                                    m_elJacobianDets, m_elIntPtCoords);
    m_elNumIntPts = (int) m_elJacobianDets.size() / m_elNum;
    gmsh::model::mesh::getBasisFunctions(m_elType[0], m_elIntType, config.elementType,
                                         m_elIntParamCoords, *new int, m_elBasisFcts);
    gmsh::model::mesh::getBasisFunctions(m_elType[0], m_elIntType, "Grad"+config.elementType,
                                         m_elIntParamCoords, *new int, m_elUGradBasisFcts);

    // Gmsh provides the derivative of the shape functions along
    // the parametric directions. We therefore compute their derivative
    // along the physical directions thanks to composed derivative.
    // The system can be expressed as J^T * df/dx = df/du
    // |dx/du dx/dv dx/dw|^T  |df/dx|   |df/du|
    // |dy/du dy/dv dy/dw|  * |df/dy| = |df/dv|
    // |dz/du dz/dv dz/dw|    |df/dz|   |df/dw|
    // NB: Instead of transposing, we take advantages of the fact
    // Lapack/Blas use column major while Gmsh provides row major.
    std::vector<double> jacobian(m_elDim*m_elDim);
    m_elGradBasisFcts.resize(m_elNum*m_elNumNodes*m_elNumIntPts*3);
    for (int el = 0; el < m_elNum; ++el) {
        for (int g = 0; g < m_elNumIntPts; ++g) {
            for (int f = 0; f < m_elNumNodes; ++f) {
                // The copy operations are not required. They're simply enforced
                // to ensure that the inputs (jacobian, grad) remains unchanged.
                for(int i=0; i<m_elDim; ++i){
                    for(int j=0; j<m_elDim; ++j) {
                        jacobian[i*m_elDim+j] = elJacobian(el, g, i, j);
                    }
                }
                std::copy(&elUGradBasisFct(g, f), &elUGradBasisFct(g, f) + m_elDim, &elGradBasisFct(el, g, f));
                eigen::solve(jacobian.data(), &elGradBasisFct(el, g, f), m_elDim);
            }
        }
    }

    assert(m_elType.size() == 1);
    assert(m_elNodeTags.size() == m_elNum*m_elNumNodes);
    assert(m_elJacobianDets.size() == m_elNum*m_elNumIntPts);
    assert(m_elBasisFcts.size() == m_elNumNodes*m_elNumIntPts);
    assert(m_elGradBasisFcts.size() == m_elNum*m_elNumIntPts*m_elNumNodes*3);

    gmsh::logger::write("------------------------------------------------");
    gmsh::logger::write("Number of Elements : " + std::to_string(m_elNum));
    gmsh::logger::write("Element dimension : " + std::to_string(m_elDim));
    gmsh::logger::write("Element Type : " + m_elName);
    gmsh::logger::write("Element Order : " + std::to_string(m_elOrder));
    gmsh::logger::write("Element Nbr Nodes : " + std::to_string(m_elNumNodes));
    gmsh::logger::write("Integration type : " + m_elIntType);
    gmsh::logger::write("Integration Nbr points : " + std::to_string(m_elNumIntPts));

    //---------------------------------------------------------------------
    // Faces
    //---------------------------------------------------------------------
    m_fDim = m_elDim - 1;
    m_fName = m_fDim == 0 ? "point"   :
              m_fDim == 1 ? "line"    :
              m_fDim == 2 ? "triangle": // Quads not yet supported.
              "None";
    m_fType = gmsh::model::mesh::getElementType(m_fName, m_elOrder);
    m_fNumNodes = m_fDim == 0 ? 1           :
                  m_fDim == 1 ? 1+m_elOrder :
                  m_fDim == 2 ? (m_elOrder+1)*(m_elOrder+2)/2 : // Triangular elements only.
                  0;

    if(m_fDim < 2)
        gmsh::model::mesh::getElementEdgeNodes(m_elType[0], m_elFNodeTags, -1);
    else
        gmsh::model::mesh::getElementFaceNodes(m_elType[0], 3, m_elFNodeTags, -1);
    m_fNumPerEl = m_elFNodeTags.size() / (m_elNum*m_fNumNodes);
    getUniqueFaceNodeTags();

    // We hereby create a single entity containing all the
    // unique faces. We call Gmsh with empty face tags and
    // retrieve directly after the auto-generated tags.
    m_fEntity = gmsh::model::addDiscreteEntity(m_fDim);
    gmsh::model::mesh::setElementsByType(m_fDim, m_fEntity, m_fType, {}, m_fNodeTags);
    m_fNodeTags.clear();
    gmsh::model::mesh::getElementsByType(m_fType, m_fTags, m_fNodeTags, m_fEntity);
    m_fNum = m_fTags.size();

    // A priori the same integration type and order is applied
    // to the surface and to the volume integrals.
    gmsh::model::mesh::getJacobians(m_fType, m_elIntType, m_fJacobians,
                                    m_fJacobianDets, m_fIntPtCoords, m_fEntity);
    m_fNumIntPts = (int) m_fJacobianDets.size() / m_fNum;
    gmsh::model::mesh::getBasisFunctions(m_fType, m_elIntType, config.elementType,
                                         m_fIntParamCoords, *new int, m_fBasisFcts);
    gmsh::model::mesh::getBasisFunctions(m_fType, m_elIntType, "Grad"+config.elementType,
                                         m_fIntParamCoords, *new int, m_fUGradBasisFcts);

    // Gmsh provides the derivative of the shape functions along
    // the parametric directions. We therefore compute their derivative
    // along the physical directions thanks to composed derivative.
    m_fGradBasisFcts.resize(m_fNum*m_fNumNodes*m_fNumIntPts*3);
    for (int f = 0; f < m_fNum; ++f) {
        for (int g = 0; g < m_fNumIntPts; ++g) {
            for (int n = 0; n < m_fNumNodes; ++n) {
                // The copy operations are not required. They're simply enforced
                // to ensure that the inputs (jacobian, grad) remains unchanged.
                for(int i=0; i<m_elDim; ++i){
                    for(int j=0; j<m_elDim; ++j) {
                        jacobian[i*m_elDim+j] = fJacobian(f, g, i, j);
                    }
                }
                std::copy(&fUGradBasisFct(g, n), &fUGradBasisFct(g, n) + m_elDim, &fGradBasisFct(f, g, n));
                eigen::solve(jacobian.data(), &fGradBasisFct(f, g, n), m_elDim);
            }
        }
    }

    // Define a normal associated to each surface.
    setFaceNormals();

    assert(m_elFNodeTags.size() == m_elNum*m_fNumPerEl*m_fNumNodes);
    assert(m_fJacobianDets.size() == m_fNum*m_fNumIntPts);
    assert(m_fBasisFcts.size() == m_fNumNodes*m_fNumIntPts);
    assert(m_fNormals.size() == m_Dim*m_fNum*m_fNumIntPts);

    gmsh::logger::write("------------------------------------------------");
    gmsh::logger::write("Number of Faces : " + std::to_string(m_fNum));
    gmsh::logger::write("Faces per Element : " + std::to_string(m_fNumPerEl));
    gmsh::logger::write("Face dimension : " + std::to_string(m_fDim));
    gmsh::logger::write("Face Type : " + m_fName);
    gmsh::logger::write("Face Nbr Nodes : " + std::to_string(m_fNumNodes));
    gmsh::logger::write("Integration Nbr points : " + std::to_string(m_fNumIntPts));

    //---------------------------------------------------------------------
    // Connectivity
    //---------------------------------------------------------------------

    // Assign corresponding faces to each element, we use the
    // fact that the node tags per face has already been ordered
    // see, 'getUniqueFaceNodeTags()'
    m_fNbrElIds.resize(m_fNum);
    for(int el=0; el<m_elNum; ++el) {
        for(int elF=0; elF<m_fNumPerEl; ++elF) {
            for(int f=0; f<m_fNum; ++f) {
                if(std::equal(&fNodeTagOrdered(f), &fNodeTagOrdered(f)+m_fNumNodes,
                              &elFNodeTagOrdered(el, elF))) {
                    m_elFIds.push_back(f);
                    m_fNbrElIds[f].push_back(el);
                }
            }
        }
    }
    assert(m_elFIds.size() == m_elNum*m_fNumPerEl);

    // For efficiency purposes we also directly store the mapping
    // between face node id and element node id. For example, the
    // 3rd node of the face correspond to the 7th of the element.
    m_fNToElNIds.resize(m_fNum);
    for(int f=0; f<m_fNum; ++f) {
        for(int nf=0; nf<m_fNumNodes; ++nf) {
            for(int el : m_fNbrElIds[f]) {
                for(int nel=0; nel<m_elNumNodes; ++nel) {
                    if(fNodeTag(f, nf) == elNodeTag(el, nel))
                        m_fNToElNIds[f].push_back(nel);
                }
            }
        }
    }

    // Define the normal direction per element with respect
    // to the normal of the face (+1 if outer,-1 if inner)
    double dotProduct;
    std::vector<double> m_elBarycenters, fNodeCoord(3), elOuterDir(3), paramCoords;
    gmsh::model::mesh::getBarycenters(m_elType[0], -1, false, true, m_elBarycenters);
    for(int el=0; el<m_elNum; ++el) {
        for(int f=0; f<m_fNumPerEl; ++f) {
            dotProduct = 0.0;
            gmsh::model::mesh::getNode(elFNodeTag(el, f), fNodeCoord, paramCoords);
            for(int x=0; x<3; x++) {
                elOuterDir[x] = fNodeCoord[x] - m_elBarycenters[el*3+x];
                dotProduct += elOuterDir[x]*fNormal(elFId(el, f), 0, x);
            }
            if(dotProduct >= 0)
                m_elFOrientation.push_back(1);
            else
                m_elFOrientation.push_back(-1);
        }
    }
    assert(m_elFOrientation.size() == m_elNum*m_fNumPerEl);

    gmsh::logger::write("------------------------------------------------");
    gmsh::logger::write("Element-Face connectivity sucessfully done.");

    //---------------------------------------------------------------------
    // Boundary conditions
    //---------------------------------------------------------------------

    // Check if a face is a boundary or not and orientate
    // the normal at boundaries in the outward direction.
    // This convention is particularly useful for BCs.
    for(int f=0; f<m_fNum; ++f){
        if(m_fNbrElIds[f].size()<2) {
            m_fIsBoundary.push_back(true);
            for(int lf=0; lf<m_fNumPerEl; ++lf) {
                if(elFId(fNbrElId(f, 0), lf) == f){
                    for(int g=0; g<m_fNumIntPts; ++g){
                        fNormal(f, g, 0) *= elFOrientation(fNbrElId(f, 0), lf);
                        fNormal(f, g, 1) *= elFOrientation(fNbrElId(f, 0), lf);
                        fNormal(f, g, 2) *= elFOrientation(fNbrElId(f, 0), lf);
                    }
                    elFOrientation(fNbrElId(f, 0), lf) = 1;
                }
            }
        }
        else {
            m_fIsBoundary.push_back(false);
        }
    }
    assert(m_fIsBoundary.size() == m_fNum);

    // Retrieve faces and nodes for boundary conditions
    m_fNeumann.resize(m_fNum);
    std::vector<int> nodeTags;
    std::vector<double> coord;
    for (auto const& physBC : config.physBCs) {
        auto physTag = physBC.first;
        auto BCtype = physBC.second.first;
        auto BCvalue = physBC.second.second;
        gmsh::model::mesh::getNodesForPhysicalGroup(m_fDim, physTag, nodeTags, coord);
        if(BCtype == "Dirichelet") {
            for(int n=0; n<nodeTags.size(); ++n) {
                for(int el=0; el<m_elNum; ++el) {
                    for(int nel=0; nel<m_elNumNodes; ++nel) {
                        if(nodeTags[n] == elNodeTag(el, nel))
                            m_elNodeDirichelet.push_back(std::make_pair(el*m_elNumNodes+nel, BCvalue));
                    }
                }
            }
        }
        else if(BCtype == "Neumann") {
            for(int f=0; f<m_fNum; ++f) {
                if(m_fIsBoundary[f] && std::find(nodeTags.begin(), nodeTags.end(), fNodeTag(f)) != nodeTags.end())
                    m_fNeumann[f] = std::make_pair(true, BCvalue);
                else
                    m_fNeumann[f] = std::make_pair(false, 0);
            }
        }
    }

    gmsh::logger::write("------------------------------------------------");
    gmsh::logger::write("Boundary conditions sucessfuly loaded.");
    gmsh::logger::write("------------------------------------------------");
}

// Precompute and store the mass matris for all elements in m_elMassMatrix
void Mesh::precomputeMassMatrix() {
    m_elMassMatrices.resize(m_elNum*m_elNumNodes*m_elNumNodes);
    for(int el=0; el<m_elNum; ++el) {
        getElMassMatrix(el, true, &elMassMatrix(el));
    }
}

// Compute the element mass matrix.
// elId : id of the element
// inverse : if the mass matrix must be inverted in place
// elMassMatrix : vector of size [NumNodes*NumNodes] containing at the
//                output the mass matrix. (row major)
void Mesh::getElMassMatrix(const int el, const bool inverse, double *elMassMatrix) {
    for(int i=0; i<m_elNumNodes; ++i) {
        for (int j=0; j<m_elNumNodes; ++j) {
            elMassMatrix[i*m_elNumNodes+j] = 0.0;
            for(int g=0; g<m_elNumIntPts; g++) {
                elMassMatrix[i*m_elNumNodes+j] += elBasisFct(g,i)*elBasisFct(g,j)*
                                                  elWeight(g)*elJacobianDet(el, g);
            }
        }
    }
    if(inverse)
        eigen::inverse(elMassMatrix, m_elNumNodes);
}

// Compute the element stiffness/convection matrix.
// elId : id of the element
// a : convection vector
// u : solution at element node
// elStiffMatrix : vector of size [NumNodes*NumNodes] containing at the
//                output the mass matrix. (row major)
void Mesh::getElStiffVector(const int el, double* a, double* u, double *elStiffVector) {
    for(int i=0; i<m_elNumNodes; ++i) {
        elStiffVector[i] = 0.0;
        for (int j = 0; j < m_elNumNodes; ++j) {
            for(int g=0; g<m_elNumIntPts; g++) {
                elStiffVector[i] += eigen::dot(a, &elGradBasisFct(el, g, i), m_Dim)*
                                    elBasisFct(g,j)*u[el*m_elNumNodes+j]*elWeight(g)*elJacobianDet(el, g);
            }
        }
    }
}

// Compute numerical flux through surface 'f'
// fId : id of the face
// a : convection vector
// u : global solution
// F : at the end return the flux at face node
void Mesh::getFlux(const int f, double* a, double * u, double* F) {
    if(m_fIsBoundary[f]) {
        if(m_fNeumann[f].first)
            std::fill(F, F+m_fNumNodes, m_fNeumann[f].second);
    }
    else {
        std::vector<double> FIntPts(m_fNumIntPts, 0);
        // Get Flux at gauss points
        for(int g=0; g<m_fNumIntPts; ++g) {
            for (int i = 0; i < m_fNumNodes; ++i) {
                // Lax-friedrich flux
                FIntPts[g] += eigen::dot(&fNormal(f, g), a, m_Dim) * fBasisFct(g, i) *
                              ((2-m_numFluxCoeff)*u[fNbrElId(f, 0) * m_elNumNodes + fNToElNId(f, i, 0)] +
                                  m_numFluxCoeff* u[fNbrElId(f, 1) * m_elNumNodes + fNToElNId(f, i, 1)]) / 2.;
            }
        }
        // Surface integral
        for(int n=0; n<m_fNumNodes; ++n) {
            F[n] = 0;
            for(int g=0; g<m_fNumIntPts; ++g) {
                F[n] += fWeight(g)*fBasisFct(g, n)*FIntPts[g]*fJacobianDet(f, g);
            }
        }
    }
}

// Precompute the flux through all surfaces
void Mesh::precomputeFlux(double* a, double * u) {
    m_fFlux.resize(m_fNum*m_fNumNodes);
    #pragma omp parallel for schedule(static) num_threads(config.numThreads)
    for(int f=0; f<m_fNum; ++f) {
        getFlux(f, a, u, &fFlux(f));
    }
}

// Compute Numerical Flux through element 'el'
void Mesh::getElFlux(const int el, double* F) {
    // Reset flux from previous iterations
    for(int n=0; n<m_elNumNodes; ++n)
        F[n] = 0;
    // Compute new Flux
    int i;
    for(int f=0; f<m_fNumPerEl; ++f) {
        el == fNbrElId(elFId(el, f), 0) ? i = 0 : i = 1;
        for(int nf=0; nf<m_fNumNodes; ++nf) {
            F[fNToElNId(elFId(el, f), nf, i)] += elFOrientation(el, f)*fFlux(elFId(el, f), nf);
        }
    }
}

// Set numerical flux type and correpsonding coefficient in lax-friedrich formula.
// Reclassify correctly neighbor elements with respect to the flux direction.
// NbrEl 0 : upstream and NbrEl 1 : downstream
void Mesh::setNumFlux(std::string fluxType, double *a, double fluxCoeff) {

    m_numFluxType = fluxType;
    if(m_numFluxType == "average")
        m_numFluxCoeff = 1;
    else if(m_numFluxType == "upwind")
        m_numFluxCoeff = 0;
    else
        m_numFluxCoeff = fluxCoeff;

    int elf;
    for(int f=0; f<m_fNum; ++f) {
        for(int lf=0; lf<m_fNumPerEl; ++lf){
            if(elFId(fNbrElId(f, 0), lf) == f)
                elf = lf;
        }
        if(!m_fIsBoundary[f]) {
            if(elFOrientation(fNbrElId(f, 0), elf)*eigen::dot(a, &fNormal(f), m_Dim) <= 0) {
                std::swap(m_fNbrElIds[f][0], m_fNbrElIds[f][1]);
                for(int nf=0; nf<m_fNumNodes; ++nf)
                    std::swap(fNToElNId(f, nf, 0), fNToElNId(f, nf, 1));
            }
        }
    }
}

void Mesh::enforceDiricheletBCs(double* u) {
    for(int n=0; n<m_elNodeDirichelet.size(); ++n)
        u[m_elNodeDirichelet[n].first] = m_elNodeDirichelet[n].second;
};

// Compute and store normal for each faces
void Mesh::setFaceNormals() {

    std::vector<double> normal(m_Dim);

    for(int f=0; f<m_fNum; ++f) {
        for(int g=0; g<m_fNumIntPts; ++g) {

            switch (m_fDim) {
                case 0: {
                    normal = {1, 0, 0};
                    break;
                }
                case 1: {
                    std::vector<double> normalPlane = {0, 0, 1};
                    eigen::cross(&fGradBasisFct(f, g, 0), normalPlane.data(), normal.data());
                    if(eigen::dot(&fGradBasisFct(f, g), &fGradBasisFct(f, 0), m_Dim) < 0) {
                        for (int x = 0; x < m_Dim; ++x) {
                            normal[x] = -normal[x];
                        }
                    }
                }
                    break;
                case 2: {
                    std::vector<double> coord1, coord2, coord3, paramCoords;
                    gmsh::model::mesh::getNode(fNodeTag(f, 0), coord1, paramCoords);
                    gmsh::model::mesh::getNode(fNodeTag(f, 1), coord2, paramCoords);
                    gmsh::model::mesh::getNode(fNodeTag(f, 2), coord3, paramCoords);
                    std::vector<double> v1 = {coord2[0]-coord1[0], coord2[1]-coord1[1], coord2[2]-coord1[2]};
                    std::vector<double> v2 = {coord3[0]-coord1[0], coord3[1]-coord1[1], coord3[2]-coord1[2]};
                    eigen::cross(v1.data(), v2.data(), normal.data());
                    break;
                }
            }
            eigen::normalize(normal.data(), m_Dim);
            m_fNormals.insert(m_fNormals.end(), normal.begin(), normal.end());
        }
    }

    // Plots
    /*std::vector<double> viewNormals;
    for(int f=0; f<m_fNum; f++) {
        for(int g=0; g<m_fNumIntPts; ++g){
            for(int x=0; x<3; ++x)
                viewNormals.push_back(fIntPtCoord(f, g, x));
            for(int x=0; x<3; ++x)
                viewNormals.push_back(fNormal(f, g, x));
        }
    }
    int normalTag = 1;
    gmsh::view::add("normals", normalTag);
    gmsh::view::addListData(normalTag, "VP", m_fNum*m_fNumIntPts, viewNormals);
    gmsh::view::write(normalTag, "normal.pos");*/

}

// Return the list of nodes for each unique face given a list of node per face and per elements
// elFNodeTags : nodes for each face of each element
//               [e1f1n1, e1f1n2, ... , e1f2n1, ..., e2f1n1, ....]
// fNodeTags : nodes for each unique face
//             [f1n1, f1n2, ... , f2n1, f2n2, ....]
void Mesh::getUniqueFaceNodeTags() {

    // Ordering per face for efficient comparison
    m_elFNodeTagsOrdered = m_elFNodeTags;

    for(int i=0; i<m_elFNodeTagsOrdered.size(); i+=m_fNumNodes)
        std::sort(m_elFNodeTagsOrdered.begin()+i, m_elFNodeTagsOrdered.begin()+(i+m_fNumNodes));

    // Unordered keep gmsh order while ordered array are used for comparison
    m_fNodeTags = m_elFNodeTags;
    m_fNodeTagsOrdered = m_elFNodeTagsOrdered;

    // Remove identical faces by comparing ordered arrays.
    std::vector<int>::iterator it_delete;
    std::vector<int>::iterator it_deleteUnordered;
    std::vector<int>::iterator it_unordered = m_fNodeTags.begin();
    for(std::vector<int>::iterator it_ordered = m_fNodeTagsOrdered.begin(); it_ordered != m_fNodeTagsOrdered.end();) {

        it_deleteUnordered = it_unordered+m_fNumNodes;
        for(it_delete=it_ordered+m_fNumNodes; it_delete != m_fNodeTagsOrdered.end(); it_delete+=m_fNumNodes) {
            if(std::equal(it_ordered, it_ordered+m_fNumNodes, it_delete))
                break;
            it_deleteUnordered+=m_fNumNodes;
        }

        if(it_delete != m_fNodeTagsOrdered.end()){
            m_fNodeTagsOrdered.erase(it_delete, it_delete+m_fNumNodes);
            m_fNodeTags.erase(it_deleteUnordered, it_deleteUnordered+m_fNumNodes);
        }
        else {
            it_ordered+=m_fNumNodes;
            it_unordered+=m_fNumNodes;
        }
    }
}