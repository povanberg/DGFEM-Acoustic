#include <string>
#include <vector>
#include <Eigen/Dense>

#ifndef DGALERKIN_FACE_H
#define DGALERKIN_FACE_H

class Face {
    private:
        //----------------------------------------------------------------------------------------------------------
        int tag;                                                            // Tag of the face
        int dim;                                                            // Dimension of the face
        int numNodes;                                                       // Number of nodes on surface
        int type;                                                           // Type of the face
        std::string name;                                                   // Face Type name
        std::vector<int> nodeTags;                                          // List of nodes of the face
        std::vector<int> elementTags;                                       // List of adjacent element tags
        Eigen::Vector3d normal;                                             // Normal (i) : n_i
        //----------------------------------------------------------------------------------------------------------
        void setNormal();                                                   // Compute normal of face

    public:
        Face(int tag, std::string name, int dim, int numNodes, int type, std::vector<int> &nodeTags);
        //----------------------------------------------------------------------------------------------------------
        const int &getTag();
        const std::vector<int> &getNodeTags();
        const Eigen::Vector3d &getNormal();
        //----------------------------------------------------------------------------------------------------------
        // Add adjacent element to the face
        Face &addElement(int tag);
};


#endif //DGALERKIN_FACE_H
