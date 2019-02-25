#include <cstdio>
#include <iostream>
#include <algorithm>
#include <gmsh.h>
#include "element.h"

int main(int argc, char **argv){

    if(argc < 2){

        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(argv[1]);

    // declares the element objects list
    std::vector<element> elementList;

    // explore the mesh: what type of 2D elements do we have?
    std::vector<int> eleTypes;
    gmsh::model::mesh::getElementTypes(eleTypes, 2);
    
    if(eleTypes.size() != 1){

        gmsh::logger::write("Hybrid meshes not handled in this example!", "error");
        return 1;
    }

    int eleType2D = eleTypes[0];
    std::string name;
    int dim, order, numNodes;
    std::vector<double> paramCoord;
    gmsh::model::mesh::getElementProperties(eleType2D, name, dim, order, numNodes, paramCoord);
    gmsh::logger::write("2D elements are of type '" + name + "' (type = " + std::to_string(eleType2D) + ") ");

    // iterate over all surfaces, get the 2D elements and create new 1D elements for all edges
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 2);

    for(std::size_t i=0; i<entities.size(); i++){

        int s = entities[i].second;
        std::vector<int> elementTags, nodeTags;

        gmsh::model::mesh::getElementsByType(eleType2D, elementTags, nodeTags,s);
        gmsh::logger::write("- " + std::to_string(elementTags.size()) + " elements in surface " + std::to_string(s));

        // get the nodes on the edges of the 2D elements
        std::vector<int> nodes;
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, nodes, s, true);

        // makes the edges uniques
        int n = nodes.size();
        int k = 0, l= 2;

        while(k < n){
            while(l < n){
                if((nodes[k]==nodes[l] && nodes[k+1]==nodes[l+1]) || (nodes[k]==nodes[l+1] && nodes[k+1]==nodes[l])){
                    nodes.erase(nodes.begin()+l);
                    nodes.erase(nodes.begin()+l);
                    n -= 2;
                }
                l += 2;
            }
            k += 2;
            l = k+2;
        }

        // create a new discrete entity of dimension 1
        int c = gmsh::model::addDiscreteEntity(1);

        // and add new 1D elements to it, for all edges
        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodes);

        // Creates a list of objects elements
        std::vector<int> tempNodes;

        for(std::size_t i=0; i<elementTags.size(); i++){
            for(std::size_t j=0; j<numNodes; j++){

                tempNodes.push_back(nodeTags[j+i*numNodes]);
            }
        
            element newElement(elementTags[i], numNodes, tempNodes);
            elementList.push_back(newElement);
            newElement.~element();
            tempNodes.clear();
        }
    }

    // iterate over all 1D elements and get integration information
    gmsh::model::mesh::getElementTypes(eleTypes, 1);
    int eleType1D = eleTypes[0];
    std::vector<double> intpts, bf;
    int numComp;
    gmsh::model::mesh::getBasisFunctions(eleType1D, "Gauss3", "IsoParametric", intpts, numComp, bf);
    gmsh::model::getEntities(entities, 1);

    for(std::size_t i=0; i<entities.size(); i++){

        int c = entities[i].second;
        std::vector<int> elementTags, nodeTags;
        gmsh::model::mesh::getElementsByType(eleType1D, elementTags, nodeTags, c);
        gmsh::logger::write("- " + std::to_string(elementTags.size()) + " elements on curve " + std::to_string(c));
        std::vector<double> jac, det, pts;
        gmsh::model::mesh::getJacobians(eleType1D, "Gauss3", jac, det, pts, c);
    }

    for(std::size_t i=0; i<elementList.size(); i++){
        std::cout << "elementList[" << i << "].m_Tag = " << elementList[i].m_Tag << std::endl;

        for(std::size_t j=0; j<elementList[i].m_numberNodes; j++){
            std::cout << "  elementList[" << i << "].m_nodes[" << j << "] = " << elementList[i].m_nodes[j] << std::endl;
        }
    }

    gmsh::finalize();
    return 0;
}
