#include <string>
#include <vector>
#include <Eigen/Dense>

#ifndef DGALERKIN_FACE_H
#define DGALERKIN_FACE_H

class Face {
    private:
        int tag;                                                            // Tag of the face
        int dim;                                                            // Dimension of the face
        int numNodes;                                                       // Number of nodes on surface
        int type;                                                           // Type of the face
        int numIntPoints;                                                   // Number of integration points
        int numBasisFcts;                                                   // Number of basis fcts
        std::string name;                                                   // Face Type name
        std::vector<int> nodeTags;                                          // List of nodes of the face
        std::vector<int> elementTags;                                       // List of adjacent element tags
        Eigen::Vector3d normal;                                             // Normal (i) : n_i
        Eigen::MatrixXd basisFcts;                                          // (f,g) : f=basis fct, g=gauss
        Eigen::MatrixXd xPoints;                                            // (g,i) : g=int point; x_i
        Eigen::MatrixXd uPoints;                                            // (g,i) : g=int point; u_i
        Eigen::VectorXd weights;                                            // (g) : g=int point
        std::vector<Eigen::Matrix3d> jacobian;                              // [g](i,j) : g=int point; dx_i/du_j
        Eigen::VectorXd detJacobian;                                        // (g) : g=int point
        std::vector<Eigen::Vector2d> elementNodes;                          // [e](n) : e=element; n=node of the face

    public:
        Face(int tag, std::string name, int dim, int numNodes, int type, std::vector<int> &nodeTags);

        bool boundary;

        void setNormal();
        void setBasisFcts();
        const int &getTag();
        const std::vector<int> &getNodeTags();
        const Eigen::Vector3d &getNormal();
        const int &getDim();

        bool hasNode(const int tag);
        void getFluxInt(Eigen::Vector2d &faceFlux, Eigen::Vector2d &uFace, const Eigen::Vector3d &a);
        int getSecondElement(const int tag1);

        // Add adjacent element to the face
        Face &addElement(int tag, std::vector<int> eNodeTags);

        // Set Jacobian
        void setJacobian(std::vector<double> &jacobian,
                         std::vector<double> &detJacobian,
                         std::vector<double> &xPoints);


};


#endif //DGALERKIN_FACE_H
