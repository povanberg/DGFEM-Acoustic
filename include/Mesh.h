#include <string>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "configParser.h"

#ifndef DGALERKIN_MESH_H
#define DGALERKIN_MESH_H

class Mesh {

public:
    Mesh(std::string name, Config config);

    // Vector access: inline optimization (c++ macro)
    // Remove function overhead, by enforcing replacement
    // -> c+11, inline is implicit inside class core
    inline int& elId(int elTag){
        return m_elIds[elTag];
    };
    inline int& elTag(int el) {
        return m_elTags[el];
    };
    // Getter : (el, n) -> node 'n' of element 'el'.
    inline int& elNodeTag(int el, int n=0) {
        return m_elNodeTags[el*m_elNumNodes+n];
    };
    // Getter : (el, g, i, j) -> dx_i/du_j(g) at int point 'g' for element 'el'.
    inline double& elJacobian(int el, int g=0, int x=0, int u=0) {
        return m_elJacobians[el*m_elNumIntPts*9 + g*9 + u*3 + x];
    };
    // Getter : (el, g) -> det of the jacobian at int point 'g' for element 'el'.
    inline double& elJacobianDet(int el, int g=0) {
        return m_elJacobianDets[el*m_elNumIntPts + g];
    };
    // Getter : (el, g, x) -> coordinate x of the int point 'g' for element 'el'.
    inline double& elIntPtCoord(int el, int g=0, int x=0) {
        return m_elIntPtCoords[el*m_elNumIntPts*3 + g*3 + x];
    };
    // Getter : (g) -> weigth associated to g-th int point.
    inline double& elWeight(int g) {
        return m_elIntParamCoords[g*4 + 3];
    };
    // Getter : (g, i) -> i-th basis fct evaluated at g-th int point.
    inline double& elBasisFct(int g, int i=0) {
        return m_elBasisFcts[g*m_elNumNodes + i];
    };
    //  Getter : (g, i, u) -> u-th derivative of the i-th basis fct evaluated at g-th int point.
    inline double& elUGradBasisFct(int g, int i=0, int u=0) {
        return m_elUGradBasisFcts[g*m_elNumNodes*3 + i*3 + u];
    };
    //  Getter : (el, g, i, x) -> x-th derivative of the i-th basis fct evaluated at g-th int point for the element 'el'.
    inline double& elGradBasisFct(int el, int g=0, int i=0, int x=0) {
        return m_elGradBasisFcts[el*m_elNumIntPts*m_elNumNodes*3 + g*m_elNumNodes*3 + i*3 + x];
    };
    // Getter : (el, f, i) -> i-th node of the f-th face for element 'el'
    inline int& elFNodeTag(int el, int f=0, int i=0) {
        return m_elFNodeTags[el*m_fNumPerEl*m_fNumNodes + f*m_fNumNodes + i];
    };
    // Getter : (f, i) -> i-th node tags of the f-th face.
    inline int& fNodeTag(int f, int i=0) {
        return m_fNodeTags[f*m_fNumNodes + i];
    };
    // Getter : (f, g, i, j) -> dx_i/du_j(g) at int point 'g' for face 'f'.
    inline double& fJacobian(int f, int g=0, int x=0, int u=0) {
        return m_fJacobians[f*m_fNumIntPts*9 + g*9 + u*3 + x];
    };
    // Getter : (f, g) -> det of the jacobian at int point 'g' for face 'f'.
    inline double& fJacobianDet(int f, int g=0) {
        return m_fJacobianDets[f*m_fNumIntPts + g];
    };
    // Getter : (f, g, x) -> coordinate x of the int point 'g' for face 'f'.
    inline double& fIntPtCoord(int f, int g=0, int x=0) {
        return m_fIntPtCoords[f*m_fNumIntPts*3 + g*3 + x];
    };
    // Getter : (g) -> weigth associated to g-th int point.
    inline double& fWeight(int g) {
        return m_fIntParamCoords[g*4 + 3];
    };
    // Getter : (g, i) -> i-th basis fct evaluated at g-th int point.
    inline double& fBasisFct(int g, int i=0) {
        return m_fBasisFcts[g*m_fNumNodes + i];
    };
    // Getter : (f, i) -> normal of face 'f'.
    inline double& fNormal(int f, int x=0) {
        return m_fNormals[f*3 + x];
    };
    // Getter : (tag) -> id
    inline int &fId(int fTag){
        return m_fIds[fTag];
    };
    // Getter : (el, f) -> f-th face tag for element 'el'
    inline int &elFId(int el, int f=0) {
        return m_elFIds[el*m_fNumPerEl +f];
    }
    // Getter (f, el) -> el-th neighbooring element of face 'f'
    inline int &fNbrElId(int f, int el=0) {
        return m_fNbrElIds[f][el];
    }
    // Getter (f, nf, nel) -> corresponding 'el' node with 'f' node.
    inline int &fNToElNId(int f, int nf=0, int el=0) {
        return m_fNToElNIds[f][nf*m_fNbrElIds[f].size() + el];
    }
    // Getter (el, f) -> orienation of the face 'f' for element 'el'
    inline int &elFOrientation( int el, int f) {
        m_elFOrientation[el*m_fNumPerEl+f];
    }
    // Getter : (el, i, j) -> i-th row and j-th column of 'el' mass matrix
    inline double &elMassMatrix(int el, int i=0, int j=0) {
        return m_elMassMatrices[el*m_elNumNodes*m_elNumNodes + i*m_elNumNodes + j];
    }
    // Getter : (f, n) -> flux through n-th of face 'f'
    inline double &fFlux(int f, int n=0) {
        return m_fFlux[f*m_fNumNodes + n];
    }
    int getNumNodes(){
        return m_elNodeTags.size();
    }
    int getElNumNodes(){
        return m_elNumNodes;
    }
    int getElNum(){
        return m_elNum;
    }

    // Compute the element mass matrix
    void getElMassMatrix(const int el, const bool inverse, double *elMassMatrix);
    // Precompute and store the mass matris for all elements in m_elMassMatrix
    void precomputeMassMatrix();
    // Compute the element stiffness/convection matrix
    void getElStiffVector(const int el, double* a, double* u, double *elStiffVector);
    // Compute Numerical Flux through surface  'f'
    void getFlux(const int f, double* a, double* u, double* F);
    // Precompute and store the flux through all surfaces
    void precomputeFlux(double* a, double * u);
    // Compute Numerical Flux through element 'el'
    void getElFlux(const int el, double* F);

private:
    std::string name;
    Config config;
    // Physics space dim
    int m_Dim = 3;
    // Dimension of the element (and the domain)
    int m_elDim;
    // Element Types (integer)
    std::vector<int> m_elType;
    // Element Type name
    std::string m_elName;
    // Element Order
    int m_elOrder;
    // Number of nodes per element
    int m_elNumNodes;
    // Number of integration points
    int m_elNumIntPts;
    // Number of elements in dim
    int m_elNum;
    // Integration type name
    std::string m_elIntType;
    // Parametric coordinates of the element
    std::vector<double> m_elParamCoord;
    // Tags of the elements
    std::vector<int> m_elTags;
    // Map the element 'tag' to the element 'ID':
    // - element tag = unique integer associated to the element by Gmsh
    // - element id = indice in 'm_el*' vector storage
    // NB: Vector instead of map guarantee very efficient O(1) lookup.
    std::vector<int> m_elIds;
    // Tags of the nodes associated to each element
    // [e1n1, e1n2, ..., e2n1, e2n2, ...]
    std::vector<int> m_elNodeTags;
    // Jacobian evaluated at each integration points : (dx/du)
    // [e1g1Jxx, e1g1Jxy, e1g1Jxz, ... e1g1Jzz, e1g2Jxx, ..., e1gGJzz, e2g1Jxx, ...]
    std::vector<double> m_elJacobians;
    // Determinants of the jacobian evaluated at each integration points
    // [e1g1DetJ, e1g2DetJ, ... e2g1DetJ, e2g2DetJ, ...]
    std::vector<double> m_elJacobianDets;
    // x, y, z coordinates of the integration points element by element.
    // [e1g1x, e1g1y, e1g1z, ... , e1gGz, e2g1x, ...]
    std::vector<double> m_elIntPtCoords;
    // u, v, w coordinates and the weight q for each integration point
    // [g1u, g1v, g1w, g1q, g2u, ...]
    std::vector<double> m_elIntParamCoords;
    // Evaluation of the basis functions at the integration points
    // [g1f1, g1f2, ..., g2f1, g2f2, ...]
    std::vector<double> m_elBasisFcts;
    // Evaluation of the derivatives of the basis functions at the integration points
    // [g1df1/du, g1df1/dv, ..., g2df1/du, g2df1/dv, ..., g1df2/du, g1df2/dv, ...]
    std::vector<double> m_elUGradBasisFcts;
    // Evaluation of the derivatives of the basis functions at the integration points
    // [e1g1df1/dx, e1g1df1/dy, ..., e1g2df1/dx, e1g2df1/dy, ..., e1g1df2/dx, e1g1df2/dy, ...]
    std::vector<double> m_elGradBasisFcts;
    // Faces ids for each element
    // [e1f1, e1f2, ..., e2f1, e2f2, ...]
    std::vector<int> m_elFIds;
    // Node tags for each face and each element
    // [e1f1n1, e1f1n2, ..., e1f2n1, e1f2n2, ..., e2f1n1, e2f1n2, ...]
    // NB: Contains duplicated faces, each element refers to its own faces.
    std::vector<int> m_elFNodeTags;
    // Contains 1 or -1, if the outward element face is in the same direction
    // as the face normal or -1 if not.
    // [e1f1, e1f2, ..., e2f1, e2f2]
    std::vector<int> m_elFOrientation;

    // Face dimension
    int m_fDim;
    // Face type name
    std::string m_fName;
    // Face type
    int m_fType;
    // Number of nodes per face
    int m_fNumNodes;
    // Number of faces per element
    int m_fNumPerEl;
    // Entity containing all the faces
    int m_fEntity;
    // Number of unique faces
    int m_fNum;
    // Number of integration points on each face
    int m_fNumIntPts;
    // Node tags for each unique face
    // [f1n1, f1n2, ..., f2n1, f2n2, ...]
    std::vector<int> m_fNodeTags;
    // Tag for each unique face
    // [f1, f2, f3, ...]
    std::vector<int> m_fTags;
    // Map the face 'tag' to the face 'ID':
    // - face tag = unique integer associated to the face by Gmsh
    // - face id = indice in 'm_f*' vector storage
    // NB: Vector instead of map guarantee very efficient O(1) lookup.
    std::vector<int> m_fIds;
    // Jacobian evaluated at each integration points : (dx/du)
    // [f1g1Jxx, f1g1Jxy, f1g1Jxz, ... f1g1Jzz, f1g2Jxx, ..., f1gGJzz, f2g1Jxx, ...]
    std::vector<double> m_fJacobians;
    // Determinants of the jacobian evaluated at each integration points
    // [f1g1DetJ, f1g2DetJ, ... f2g1DetJ, f2g2DetJ, ...]
    std::vector<double> m_fJacobianDets;
    // x, y, z coordinates of the integration points for each faces
    // [f1g1x, f1g1y, f1g1z, ... , f1gGz, f2g1x, ...]
    std::vector<double> m_fIntPtCoords;
    // u, v, w coordinates and the weight q for each integration point
    // [g1u, g1v, g1w, g1q, g2u, ...]
    std::vector<double> m_fIntParamCoords;
    // Evaluation of the basis functions at the integration points
    // [g1f1, g1f2, ..., g2f1, g2f2, ...]
    std::vector<double> m_fBasisFcts;
    // Normal for each surface
    // [f1Nx, f1Ny, f1Nz, f2Nx, ...]
    std::vector<double> m_fNormals;
    // Is Face a boundary
    std::vector<bool> m_fIsBoundary;
    // Id of element of each side of the face
    // NB: As this number is variable, for example, the face at the boundary
    //     have one less neighbor, it requested the use of a 2D vector (or map)
    std::vector<std::vector<int>> m_fNbrElIds;
    // Map face node Ids to element node Ids
    // [f1n1e1n, f1n1e2n, ..., f1n2e1n, f1n2e2n, ....]
    // [f2n1e1n, f2n1e2n, ..., f2n2e1n, f2n2e2n, ....]
    // [                  ...                        ]
    // [fNn1e1n, fNn1e2n, ..., fNn2e1n, fNn2e2n, ....]
    std::vector<std::vector<int>> m_fNToElNIds;

    // Element mass matrix stored contiguously (row major)
    // [e1m11, e1m12, ..., e1m21, e1m22, ..., e2m11, ...]
    std::vector<double> m_elMassMatrices;
    // Flux through all faces
    // [f1n1, f1n2, ..., f2n1, f2n2, ...]
    std::vector<double> m_fFlux;
};

#endif //DGALERKIN_MESH_H