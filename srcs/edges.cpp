#include <cstdio>
#include <iostream>
#include <algorithm>
#include <gmsh.h>

int main(int argc, char **argv){

    if(argc < 2){
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(argv[1]);

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

        // sorts the edges and makes them uniques
        int n = nodes.size();
        int j = 0, k= 2;

        while(j < n){
            while(k < n){
                if((nodes[j]==nodes[k] && nodes[j+1]==nodes[k+1]) || (nodes[j]==nodes[k+1] && nodes[j+1]==nodes[k])){
                    nodes.erase(nodes.begin()+k);
                    nodes.erase(nodes.begin()+k);
                    n -= 2;
                }
                k += 2;
            }
            j += 2;
            k = j+2;
        }

        // create a new discrete entity of dimension 1
        int c = gmsh::model::addDiscreteEntity(1);

        // and add new 1D elements to it, for all edges
        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodes);

        // this could be enriched with additional info: each topological edge could
        // be associated with the tag of its parent element; in the sorting process
        // eliminating duplicates a second tag can be associated for internal edges,
        // allowing to keep track of neighbors
    }

    //gmsh::write("edges.msh");

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

    //gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
