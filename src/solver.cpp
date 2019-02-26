#include "solver.h"
#include <configParser.h>
#include <vector>
#include <Mesh.h>
#import "logger.h"

void zeroInitializer(std::vector<std::vector<double>> &solution, const Mesh &mesh) {
    for (unsigned int i=0; i<mesh.faceNodeTags.size(); ++i) {
        std::vector<double> initValue(1,0);
        solution.push_back(initValue);
    }
}

double integrateGauss(){
    return 0.;
}

namespace solver {

    void solveForwardEuler(Mesh &mesh, const Config config) {

        // 0: Current solution
        // 1: Next iteration solution
        std::vector<std::vector<double>> solution;

        // 1: Set Initial conditions
        zeroInitializer(solution, mesh);

        // 2: Precompute
        //auto elementTag = mesh.elementTags[0];
        //error("element: %i", elementTag);

        // 2.a: Mass matrix: [(row * columns) + column]
        //std::vector<double> massMatrix(mesh.elementNumNodes*mesh.elementNumNodes);

        /*std::vector<double> integrationPoints;
        std::vector<double> basisFunctions;
        int numComponents;
        gmsh::model::mesh::getBasisFunctions(
                mesh.elementType,
                "Gauss3",
                "Lagrange",
                integrationPoints,
                numComponents,
                basisFunctions);

        int numIntPoints = (int) integrationPoints.size() / 4;
        // For each basis function combination
        for(unsigned int i=0; i<mesh.elementNumNodes; ++i) {
            for (unsigned int j = 0; j < mesh.elementNumNodes; ++j) {

            }
        }*/







        // Save results example
        /*int viewTag = gmsh::view::add("nameConfig");
        double t0 = 0;
        double tf = 1;
        int step = 0;
        double timeStep = 0.1;
        for(double t=t0; t<tf; t+=timeStep, ++step) {
            for(unsigned int i=0; i<solution.size();++i)
            {
                solution[i][0] += 10;
            }
            log("[%f/%f] Simulation time step.", t, tf);

            std::vector<std::string> names;
            gmsh::model::list(names);

            std::vector<int> nodeTags;
            std::vector<double> coords;
            std::vector<double> parametricCoord;
            gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoord);

            gmsh::view::addModelData(viewTag, step, names[0], "NodeData", mesh.faceNodeTags, solution, t, 1);
        }

        gmsh::view::write(viewTag, "data.msh");*/
    }
}