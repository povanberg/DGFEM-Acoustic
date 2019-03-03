#include "solver.h"
#include <configParser.h>
#include <vector>
#include <Mesh.h>
#include "logger.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

double integrateGauss(){
    return 0.;
}

namespace solver {


    void solveForwardEuler(Mesh &mesh, const Config config) {

        // 0: Current solution
        // 1: Next iteration solution
        std::vector<std::vector<double>> solution;

        Eigen::SparseMatrix<double> M;
        mesh.getMassMatrix(M);

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