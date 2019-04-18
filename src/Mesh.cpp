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

    // Reclassify neighbor element, the first one has the
    // normal oriented in the same direction as the face.
    int elf;
    for(int f=0; f<m_fNum; ++f) {
        for(int lf=0; lf<m_fNumPerEl; ++lf){
            if(elFId(fNbrElId(f, 0), lf) == f)
                elf = lf;
        }
        if(m_fNbrElIds.size() == 2) {
            if(elFOrientation(fNbrElId(f, 0), elf) <= 0) {
                std::swap(m_fNbrElIds[f][0], m_fNbrElIds[f][1]);
                for(int nf=0; nf<m_fNumNodes; ++nf)
                    std::swap(fNToElNId(f, nf, 0), fNToElNId(f, nf, 1));
            }
        }
    }

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
    m_fBC.resize(m_fNum);
    std::vector<int> nodeTags;
    std::vector<double> coord;
    for (auto const& physBC : config.physBCs) {
        auto physTag = physBC.first;
        auto BCtype = physBC.second.first;
        auto BCvalue = physBC.second.second;
        gmsh::model::mesh::getNodesForPhysicalGroup(m_fDim, physTag, nodeTags, coord);
        if(BCtype == "Reflecting") {
            for(int f=0; f<m_fNum; ++f) {
                if(m_fIsBoundary[f] && std::find(nodeTags.begin(), nodeTags.end(), fNodeTag(f)) != nodeTags.end())
                    m_fBC[f] = 1;
            }
        }
        else {
            for(int f=0; f<m_fNum; ++f) {
                if(m_fIsBoundary[f] && std::find(nodeTags.begin(), nodeTags.end(), fNodeTag(f)) != nodeTags.end())
                    m_fBC[f] = 0;
            }
        }
    }

    gmsh::logger::write("------------------------------------------------");
    gmsh::logger::write("Boundary conditions sucessfuly loaded.");
    gmsh::logger::write("------------------------------------------------");

    m_fFlux.resize(m_fNum*m_fNumNodes);
    uGhost = std::vector<std::vector<double>>(4,
             std::vector<double>(getNumNodes()));
    FluxGhost = std::vector<std::vector<std::vector<double>>>(4,
                std::vector<std::vector<double>>(getNumNodes(),
                std::vector<double>(3)));
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
void Mesh::getElStiffVector(const int el, std::vector<std::vector<double>> &Flux, std::vector<double> &u, double *elStiffVector) {
    for(int i=0; i<m_elNumNodes; ++i) {
        elStiffVector[i] = 0.0;
        for (int j = 0; j < m_elNumNodes; ++j) {
            for(int g=0; g<m_elNumIntPts; g++) {
                elStiffVector[i] += eigen::dot(Flux[el*m_elNumNodes+j].data(), &elGradBasisFct(el, g, i), m_Dim)*
                                    elBasisFct(g,j)*elWeight(g)*elJacobianDet(el, g);
            }
        }
    }
}

// Precompute the flux through all surfaces
void Mesh::precomputeFlux(std::vector<double> &u, std::vector<std::vector<double>> &Flux, int eq) {

    #pragma omp parallel num_threads(config.numThreads)
    {
        int elUp, elDn;
        std::vector<double> FIntPts(m_fNumIntPts, 0);
        std::vector<double> Fnum(m_Dim, 0);
        #pragma omp parallel for schedule(static)
        for(int f=0; f<m_fNum; ++f) {
            std::fill(FIntPts.begin(), FIntPts.end(), 0);
            if (m_fIsBoundary[f]) {
                if (m_fBC[f] == 1) {
                    for (int i = 0; i < m_fNumNodes; ++i) {
                        elUp = fNbrElId(f, 0) * m_elNumNodes + fNToElNId(f, i, 0);
                        for (int g = 0; g < m_fNumIntPts; ++g) {
                            for (int x = 0; x < m_Dim; ++x)
                                Fnum[x] = (Flux[elUp][x] + FluxGhost[eq][elUp][x]) / 2.;
                            FIntPts[g] += eigen::dot(&fNormal(f, g), Fnum.data(), m_Dim) * fBasisFct(g, i);
                        }
                    }
                } else {
                    for (int i = 0; i < m_fNumNodes; ++i) {
                        elUp = fNbrElId(f, 0) * m_elNumNodes + fNToElNId(f, i, 0);
                        for (int g = 0; g < m_fNumIntPts; ++g) {
                            FIntPts[g] += eigen::dot(&fNormal(f, g), &Flux[elUp][0], m_Dim) * fBasisFct(g, i);
                        }
                    }
                }
            } else {
                for (int i = 0; i < m_fNumNodes; ++i) {
                    elUp = fNbrElId(f, 0) * m_elNumNodes + fNToElNId(f, i, 0);
                    elDn = fNbrElId(f, 1) * m_elNumNodes + fNToElNId(f, i, 1);
                    for (int g = 0; g < m_fNumIntPts; ++g) {
                        for (int x = 0; x < m_Dim; ++x)
                            Fnum[x] = ((Flux[elUp][x] + Flux[elDn][x]) +
                                       config.c0 * fNormal(f, g, x) * (u[elDn] - u[elUp])) / 2.;
                        FIntPts[g] += eigen::dot(&fNormal(f, g), Fnum.data(), m_Dim) * fBasisFct(g, i);
                    }
                }
            }
            // Surface integral
            for (int n = 0; n < m_fNumNodes; ++n) {
                fFlux(f, n) = 0;
                for (int g = 0; g < m_fNumIntPts; ++g) {
                    fFlux(f, n) += fWeight(g) * fBasisFct(g, n) * FIntPts[g] * fJacobianDet(f, g);
                }
            }
        }
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

void Mesh::updateFlux(std::vector<std::vector<double>> &u, std::vector<std::vector<std::vector<double>>> &Flux,
                      std::vector<double> &v0, double c0, double rho0) {

    int i, fId;
    for(int el=0; el<m_elNum; ++el){
        for(int n=0; n<m_elNumNodes; ++n) {
            i = el*m_elNumNodes + n;

            // Pressure flux
            Flux[0][i] = {v0[0]*u[0][i] + rho0*c0*c0*u[1][i],
                          v0[1]*u[0][i] + rho0*c0*c0*u[2][i],
                          v0[2]*u[0][i] + rho0*c0*c0*u[3][i]};
            // Vx
            Flux[1][i] = {v0[0]*u[1][i] + u[0][i]/rho0,
                          v0[1]*u[1][i],
                          v0[2]*u[1][i]};
            // Vy
            Flux[2][i] = {v0[0]*u[2][i],
                          v0[1]*u[2][i] + u[0][i]/rho0,
                          v0[2]*u[2][i]};
            // Vz
            Flux[3][i] = {v0[0]*u[3][i],
                          v0[1]*u[3][i],
                          v0[2]*u[3][i] + u[0][i]/rho0};

        }

        // Ghost elements
        for(int f=0; f<m_fNumPerEl; ++f) {
            fId = elFId(el, f);
            if(m_fIsBoundary[fId]) {
                for(int n=0; n<m_elNumNodes; ++n) {
                    i = el*m_elNumNodes + n;
                    double dot = fNormal(fId,0,0)*u[1][i] +
                                 fNormal(fId,0,1)*u[2][i] +
                                 fNormal(fId,0,2)*u[3][i];

                    uGhost[0][i] = u[0][i];
                    uGhost[1][i] = u[1][i] - 2*dot*fNormal(fId,0,0);
                    uGhost[2][i] = u[2][i] - 2*dot*fNormal(fId,0,1);
                    uGhost[3][i] = u[3][i] - 2*dot*fNormal(fId,0,2);

                    // Pressure flux
                    FluxGhost[0][i] = {v0[0]*uGhost[0][i] + rho0*c0*c0*uGhost[1][i],
                                       v0[1]*uGhost[0][i] + rho0*c0*c0*uGhost[2][i],
                                       v0[2]*uGhost[0][i] + rho0*c0*c0*uGhost[3][i]};
                    // Vx
                    FluxGhost[1][i] = {v0[0]*uGhost[1][i] + uGhost[0][i]/rho0,
                                       v0[1]*uGhost[1][i],
                                       v0[2]*uGhost[1][i]};
                    // Vy
                    FluxGhost[2][i] = {v0[0]*uGhost[2][i],
                                       v0[1]*uGhost[2][i] + uGhost[0][i]/rho0,
                                       v0[2]*uGhost[2][i]};
                    // Vz
                    FluxGhost[3][i] = {v0[0]*uGhost[3][i],
                                       v0[1]*uGhost[3][i],
                                       v0[2]*uGhost[3][i] + uGhost[0][i]/rho0};
                }

            }
        }
    }

}

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
                    break;
                }
                case 2: {
                    eigen::cross(&fGradBasisFct(f, g, 0), &fGradBasisFct(f, g, 1), normal.data());
                    break;
                }
            }
            eigen::normalize(normal.data(), m_Dim);
            m_fNormals.insert(m_fNormals.end(), normal.begin(), normal.end());
        }
    }
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
