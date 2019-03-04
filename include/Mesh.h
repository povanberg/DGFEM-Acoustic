#include <gmsh.h>
#include <map>
#include "Element.h"
#include <Eigen/Sparse>

#ifndef DGALERKIN_MESH_H
#define DGALERKIN_MESH_H

class Mesh {

    private:

        int dim;                                                    // Mesh dimension
        std::string name;                                           // File name
        std::map<std::string, std::pair<int,int>> physDimTags;      // Associate physical name to tag/dim
        std::vector<int> domainEntityTags;                          // List of all geometrical entities in domain
        std::vector<int> diricheletEntityTags;                      // List of all geometrical entities in dirichelet

    public:

        Mesh(std:: string name);

        std::vector<Element> elements;                              // List of elements in Mesh
        Element &getElement(const int tag);

        void getMassMatrix(Eigen::SparseMatrix<double> &M);
        void getStiffMatrix(Eigen::SparseMatrix<double> &K, const Eigen::Vector3d &a);
        void getFlux(Eigen::VectorXd &F, const Eigen::Vector3d &a);

        Element &getDataElement(const int tag);
        void getNodeTags(std::vector<int> &nodeTags);
        void getData(Eigen::VectorXd &data);
        void setData(Eigen::VectorXd &data);

};


#endif //DGALERKIN_MESH_H
