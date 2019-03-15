#include <string>
#include <gmsh.h>
#include <assert.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "configParser.h"
#include "Mesh.h"
#include "utils.h"

// Return a unique list nodes for each face given a list
// of faces per elements.
void getUniqueFaceNodeTags(std::vector<int> &elFNodeTags, const int fNumPerEl, std::vector<int> &fNodeTags) {

    // Sort nodes on each faces
    for(int i=0; i<elFNodeTags.size(); i+=fNumPerEl)
        std::sort(elFNodeTags.begin()+i, elFNodeTags.begin()+(i+fNumPerEl));

    // Using assignment operator to copy one vector to other
    fNodeTags = elFNodeTags;

    // Remove identical faces
    for(std::vector<int>::iterator it = fNodeTags.begin(); it != fNodeTags.end();) {

        std::vector<int>::iterator it_d;
        for(it_d= it+fNumPerEl; it_d != fNodeTags.end(); it_d+=fNumPerEl) {
            if(std::equal(it, it+fNumPerEl, it_d))
                break;
        }

        if(it_d != fNodeTags.end())
            fNodeTags.erase(it_d, it_d+fNumPerEl);
        else
            it+=fNumPerEl;
    }
}

#include <iomanip>
template<typename Container>
void print(const Container& cont, int row = 1, bool colMajor=false) {

    if(colMajor){
        for(int rowIt=0; rowIt<row; ++rowIt){
            int colIt = 0;
            for (auto const& x : cont) {
                if(colIt%row == rowIt) {
                    std::cout << std::setprecision(4) << std::left << std::setw(10) << x << " ";
                }
                colIt++;
            }
            std::cout << std::endl;
        }
    }
    else {
        int colIt = 0;
        for (auto const& x : cont) {
            std::cout << std::setprecision(4) << std::left << std::setw(10) << x << " ";
            colIt++;
            if(colIt%row == 0)
                std::cout << std::endl;
        }
    }

}

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
    // polynomials of order n with p=2n-1 integration points.
    m_elNum = m_elTags.size();
    m_elIntType = "Gauss" + std::to_string(m_elOrder+3);
    gmsh::model::mesh::getJacobians(m_elType[0], m_elIntType, m_elJacobians,
                                    m_elJacobianDets, m_elIntPtCoords);
    m_elNumIntPts = (int) m_elJacobianDets.size() / m_elNum;
    gmsh::model::mesh::getBasisFunctions(m_elType[0], m_elIntType, config.elementType,
                                         m_elIntParamCoords, *new int, m_elBasisFcts);
    gmsh::model::mesh::getBasisFunctions(m_elType[0], m_elIntType, "Grad"+config.elementType,
                                         m_elIntParamCoords, *new int, m_elUGradBasisFcts);

    // Gmsh provides the derivative of the shape functions alongnumNodes
    // the parametric directions. We therefore compute their derivative
    // along the physical directions thanks to composed derivative.
    // The system can be expressed as J^T * df/dx = df/du
    // |dx/du dx/dv dx/dw|^T  |df/dx|   |df/du|
    // |dy/du dy/dv dy/dw|  * |df/dy| = |df/dv|
    // |dw/du dw/dv dw/dw|    |df/dz|   |df/dw|
    // NB: Instead of transposing, we take advantages of the fact
    // Lapack/Blas use column major while Gmsh provides row major.
    int jacobianSize = 3;
    std::vector<double> jacobian(jacobianSize*jacobianSize);
    m_elGradBasisFcts = std::vector<double>(m_elNum*m_elNumNodes*m_elNumIntPts*3);
    for (int el = 0; el < m_elNum; ++el) {
        for (int g = 0; g < m_elNumIntPts; ++g) {
            for (int f = 0; f < m_elNumNodes; ++f) {
                // The copy operations are not required. They're simply enforced
                // to ensure that the inputs (jacobian, grad) remains unchanged.
                std::copy(&elJacobian(el, g), &elJacobian(el, g) + 9, jacobian.begin());
                std::copy(&elUGradBasisFct(g, f), &elUGradBasisFct(g, f) + 3, &elGradBasisFct(el, g, f));
                lapack::solve(jacobian.data(), &elGradBasisFct(el, g, f), jacobianSize);
            }
        }
    }

    // Map the element 'tag' to the face 'ID':
    // - element tag = unique integer associated to the element by Gmsh
    // - element id = indice in 'm_f*' vector storage
    // This hypothesis rely on fact that Gmsh could provide non-zero
    // starting or non-continuous sequence of face tags.
    // e.g. m_elTags = [23, 24, 25, 78, 79 ...] => m_elIds = [0, 1, 2, 3, 4]
    for(int elId=0; elId<m_elNum; ++elId) {
        if(m_elIds.size() <= m_elTags[elId])
            m_elIds.resize(m_elTags[elId]);
        m_elIds.insert(m_elIds.begin()+m_elTags[elId],elId);
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
                  m_fDim == 2 ? 3*m_elOrder : // Triangular elements only and up to order 2.
                  0;

    if(m_fDim < 2)
        gmsh::model::mesh::getElementEdgeNodes(m_elType[0], m_elFNodeTags, -1);
    else
        gmsh::model::mesh::getElementFaceNodes(m_elType[0], 3, m_elFNodeTags, -1);
    m_fNumPerEl = m_elFNodeTags.size() / (m_elNum*m_fNumNodes);
    getUniqueFaceNodeTags(m_elFNodeTags, m_fNumNodes, m_fNodeTags);


    // We hereby create a single entity containing all the
    // unique faces. We call Gmsh with empty face tags and
    // retrieve directly after the auto-generated tags.
    m_fEntity = gmsh::model::addDiscreteEntity(m_fDim);
    gmsh::model::mesh::setElementsByType(m_fDim, m_fEntity, m_fType, {}, m_fNodeTags);
    m_fNodeTags.clear();
    gmsh::model::mesh::getElementsByType(m_fType, m_fTags, m_fNodeTags, m_fEntity);
    m_fNum = m_fTags.size();

    // Map the face 'tag' to the face 'ID':
    // - face tag = unique integer associated to the face by Gmsh
    // - face id = indice in 'm_f*' vector storage
    // This hypothesis rely on fact that Gmsh could provide non-zero
    // starting or non-continuous sequence of face tags.
    // e.g. m_fTags = [23, 24, 25, 78, 79 ...] => m_fIds = [0, 1, 2, 3, 4]
    for(int fId=0; fId<m_fNum; ++fId) {
        if(m_fIds.size() <= m_fTags[fId])
            m_fIds.resize(m_fTags[fId]);
        m_fIds.insert(m_fIds.begin()+m_fTags[fId],fId);
    }

    // A priori the same integration type and order is applied
    // to the surface and to the volume integrals.
    gmsh::model::mesh::getJacobians(m_fType, m_elIntType, m_fJacobians,
                                    m_fJacobianDets, m_fIntPtCoords, m_fEntity);
    m_fNumIntPts = (int) m_fJacobianDets.size() / m_fNum;
    gmsh::model::mesh::getBasisFunctions(m_fType, m_elIntType, config.elementType,
                                         m_fIntParamCoords, *new int, m_fBasisFcts);

    // Define a normal associated to each surface. For now
    // the faces are assumed to be straight.
    int size = 3;
    std::vector<double> normal, coord1, coord2, paramCoords;
    for(int f=0; f<m_fNum; ++f) {
        switch(m_fDim) {
            case 0:
                normal = {1, 0, 0};
                break;
            case 1:
                gmsh::model::mesh::getNode(fNodeTag(f, 0), coord1, paramCoords);
                gmsh::model::mesh::getNode(fNodeTag(f, 1), coord2, paramCoords);
                normal = {coord1[1]-coord2[1], coord2[0]-coord1[0], 0};
                break;
            case 2:
                // Todo(19/03/2019)
                break;
        }
        lapack::normalize(normal.data(), size);
        m_fNormals.insert(m_fNormals.end(), normal.begin(), normal.end());
    }

    assert(m_elFNodeTags.size() == m_elNum*m_fNumPerEl*m_fNumNodes);
    assert(m_fJacobianDets.size() == m_fNum*m_fNumIntPts);
    assert(m_fBasisFcts.size() == m_fNumNodes*m_fNumIntPts);
    assert(m_fNormals.size() == m_fNum*m_fNumPerEl);

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
                if(std::equal(&fNodeTag(f), &fNodeTag(f)+m_fNumNodes,
                              &elFNodeTag(el, elF))) {
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
    std::vector<double> m_elBarycenters, fNodeCoord(3), elOuterDir(3);
    gmsh::model::mesh::getBarycenters(m_elType[0], -1, false, true, m_elBarycenters);
    for(int el=0; el<m_elNum; ++el) {
        for(int f=0; f<m_fNumPerEl; ++f) {
            dotProduct = 0.0;
            gmsh::model::mesh::getNode(elFNodeTag(el, f), fNodeCoord, paramCoords);
            for(int x=0; x<3; x++) {
                elOuterDir[x] = fNodeCoord[x] - m_elBarycenters[el*3+x];
                dotProduct += elOuterDir[x]*fNormal(elFId(el, f), x);
            }
            if(dotProduct >= 0)
                m_elFOrientation.push_back(1);
            else
                m_elFOrientation.push_back(-1);
        }
    }
    assert(m_elFOrientation.size() == m_elNum*m_fNumPerEl);

    // Check if a face is a boundary or not.
    for(int f=0; f<m_fNum; ++f){
        if(m_fNbrElIds[f].size()<2)
            m_fIsBoundary.push_back(true);
        else
            m_fIsBoundary.push_back(false);
    }
    assert(m_fIsBoundary.size() == m_fNum);


        //auto it_duplicate = std::search(it+fNumPerEl, fNodeTags.end(), it, it+fNumPerEl);
    /*print(m_fTags, 1, true);
    std::cout << "------------------------" << std::endl;
    print(m_fNodeTags, m_fNumNodes, true);
    std::cout << "------------------------" << std::endl;
    print(m_elFNodeTags, m_fNumNodes*m_fNumPerEl, true);
    std::cout << "------------------------" << std::endl;
    print(m_elFIds, m_fNumPerEl, true);*/

    /*std::vector<double> u(m_elNum*m_elNumNodes,1);
    std::vector<double> a = {1, 0, 0};
    std::vector<double> F(3);
    precomputeMassMatrix();
    precomputeFlux(a.data(), u.data());
    print(m_elNodeTags, m_elNumNodes, true);
    for(int el=0; el<m_elNum; ++el){
        getElFlux(el, F.data());
        print(F, 1, true);
        std::cout << "------------------------" << std::endl;
    }*/

}

// Precompute and store the mass matris for all elements in m_elMassMatrix
void Mesh::precomputeMassMatrix() {
    m_elMassMatrices.resize(m_elNum*m_elNumNodes*m_elNumNodes);
    for(int el=0; el<m_elNum; ++el)
        getElMassMatrix(el, true, &elMassMatrix(el));
}

// Compute the element mass matrix.
// elId : id of the element
// inverse : if the mass matrix must be inverted in place
// elMassMatrix : vector of size [NumNodes*NumNodes] containing at the
//                output the mass matrix. (row major)
void Mesh::getElMassMatrix(const int el, const bool inverse, double *elMassMatrix) {
    for(int i=0; i<m_elNumNodes; ++i) {
        for (int j = 0; j < m_elNumNodes; ++j) {
            elMassMatrix[i*m_elNumNodes+j] = 0.0;
            for(int g=0; g<m_elNumIntPts; g++) {
                elMassMatrix[i*m_elNumNodes+j] += elBasisFct(g,i)*elBasisFct(g,j)*
                                                  elWeight(g)*elJacobianDet(el, g);
            }
        }
    }
    if(inverse)
        lapack::inverse(elMassMatrix, m_elNumNodes);
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
                elStiffVector[i] += lapack::dot(a, &elGradBasisFct(el, g, i), m_Dim)*
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
    if(m_fNbrElIds[f].size()<2) {
        for(int n=0; n<m_fNumNodes; ++n)
            F[n] = 0;
    }
    else {
        // Some precomputation
        double dot = lapack::dot(&fNormal(f), a, m_Dim);
        std::vector<double> FIntPts(m_fNumIntPts, 0);
        // Surface integral
        for(int n=0; n<m_fNumNodes; ++n) {
            F[n] = 0;
            std::fill(FIntPts.begin(), FIntPts.end(), 0);
            for(int g=0; g<m_fNumIntPts; ++g) {
                for(int i=0; i<m_fNumNodes; i++)
                    FIntPts[g] += dot*fBasisFct(g, i)*(u[fNbrElId(f, 0)*m_elNumNodes+fNToElNId(f, i, 0)]+
                                                       u[fNbrElId(f, 1)*m_elNumNodes+fNToElNId(f, i, 1)])/2.;
                F[n] += fWeight(g)*fBasisFct(g, n)*FIntPts[g]*fJacobianDet(f, g);
            }
        }
    }
}

// Precompute the flux through all surfaces
void Mesh::precomputeFlux(double* a, double * u) {
    m_fFlux.resize(m_fNum*m_fNumNodes);
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


// Small code to vizualize the normals.
/*std::vector<double> viewNormals;
for(int f=0; f<m_fNum; f++) {
std::vector<double> fNodeCoord, fNodeCoord2, paramCoords;
gmsh::model::mesh::getNode(fNodeTag(f), fNodeCoord, paramCoords);
gmsh::model::mesh::getNode(fNodeTag(f, 1), fNodeCoord2, paramCoords);
for(int x=0; x<3; ++x)
viewNormals.push_back((fNodeCoord[x]+fNodeCoord2[x])/2);
for(int x=0; x<3; ++x)
viewNormals.push_back(fNormal(f, x));
}
int normalTag;
gmsh::view::add("normals", normalTag);
gmsh::view::addListData(normalTag, "VP", m_elNum, viewNormals);
gmsh::view::write(normalTag, "normal.pos");*/