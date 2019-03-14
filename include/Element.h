#include <string>
#include <vector>
#include "Face.h"
#include <Eigen/Dense>

#ifndef DGALERKIN_ELEMENT_H
#define DGALERKIN_ELEMENT_H


class Element {
    private:
        //----------------------------------------------------------------------------------------------------------
        int dim;                                                        // Element dimension
        int tag;                                                        // Element tag
        int type;                                                       // Element type
        int order;                                                      // Element order
        int numNodes;                                                   // Num nodes per Elements
        std::string name;                                               // Name of element type
        std::vector<double> paramCoord;                                 // Parametric coordinates
        std::vector<int> nodeTags;                                      // List of node tags
        //----------------------------------------------------------------------------------------------------------
        int numIntPoints;                                               // Number of Integration points
        int numBasisFcts;                                               // Number of Basis functions
        std::vector<Eigen::Matrix3d> jacobian;                          // [g](i,j) : g=int point; dx_i/du_j
        std::vector<Eigen::Matrix3d> invJacobian;                       // [g](i,j) : g=int point; du_i/dx_j
        Eigen::VectorXd detJacobian;                                    // (g) : g=int point
        Eigen::MatrixXd xPoints;                                        // (g,i) : g=int point; x_i
        Eigen::MatrixXd uPoints;                                        // (g,i) : g=int point; u_i
        Eigen::VectorXd weights;                                        // (g) : g=int point
        std::vector<Eigen::VectorXd> basisFcts;                         // [f](g) : f= basis fct; g=int point
        std::vector<Eigen::MatrixXd> gradBasisFcts;                     // [f](g,i) : f= fct; g=int point; i=df/dx_i
        std::vector<Eigen::Vector3d> normals;
        Eigen::Vector3d barycenter;
        //----------------------------------------------------------------------------------------------------------
        Eigen::VectorXd u;
public:
        Element(int dim, int tag, std::vector<int> &nodeTags, Eigen::Vector3d &barycenter);
        std::vector<Face> faces;
        //----------------------------------------------------------------------------------------------------------
        const int &getTag();
        const int &getType();
        void setNormals();
        const int &getOrder();
        const int &getNumNodes();
        std::vector<int> &getNodeTags();
        const std::string &getName();
        const Eigen::Matrix3d &getJacobian(const int g);
        const Eigen::Matrix3d &getInvJacobian(const int g);
        const Eigen::VectorXd &getDetJacobian();
        const Eigen::MatrixXd &getXPoints();
        const Eigen::MatrixXd &getUPoints();
        const Eigen::VectorXd &getWeigths();
        const Eigen::VectorXd &getBasisFcts(const int f);
        const Eigen::MatrixXd &getGradBasisFcts(const int f);
        void getMassMatrix(Eigen::MatrixXd &massMAtrix);
        void getStiffMatrix(Eigen::MatrixXd &stiffMatrix, const Eigen::Vector3d &a);
        void getData(Eigen::VectorXd &data);
        void setData(Eigen::VectorXd &data);
        //----------------------------------------------------------------------------------------------------------

        bool hasNode(const int tag);
        void getFlux(Eigen::VectorXd &Flux, const Eigen::Vector3d &a, std::vector<Element> &elements, Eigen::VectorXd &u);

        // Add face to element
        Element &addFace(Face face);

        // Set jacobian for the element
        Element &setJacobian(std::vector<double> &jacobian,
                             std::vector<double> &detJacobian,
                             std::vector<double> &xPoints,
                             int numIntPoints);

        // Set basis functions for the element
        Element &setBasis(std::vector<double> &ubasisFct,
                          std::vector<double> &ugradBasisFct,
                          std::vector<double> &uPoints);
};


#endif //DGALERKIN_ELEMENT_H
