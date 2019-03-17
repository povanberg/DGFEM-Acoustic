#include <vector>
#include <gmsh.h>
#include <utils.h>
#include <iostream>
#include "configParser.h"
#include "Mesh.h"

namespace solver {

    // LAPack
    extern "C" void dgemv_(char& TRANS, int& M, int& N, double& a, double* A,
                           int& LDA, double* X, int& INCX, double& beta, double* Y, int& INCY);

    // Initial confiditon
    double f(std::vector<double> x) {
        //return (1) * exp(-((x[0] + 3.2) * (x[0] + 3.2))/(0.5));
        return (1) * exp(-((x[0] - 0) * (x[0] - 0) + (x[1]- 0) * (x[1]- 0))/0.5);
    }

    void solveForwardEuler(Mesh &mesh, Config config) {

        // Solution vectors
        int N = mesh.getElNumNodes();
        std::vector<int> nodeTags(&mesh.elNodeTag(0), &mesh.elNodeTag(0) + mesh.getNumNodes());
        std::vector<double> u(mesh.getNumNodes());
        std::vector<double> u_prev(mesh.getNumNodes());

        // Save results initialization
        int viewTag = gmsh::view::add("Results");
        std::vector<std::string> names;
        gmsh::model::list(names);
        std::vector<std::vector<double>> solution(nodeTags.size(), std::vector<double>(1));

        // Set Initial conditions
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        for(int n=0; n<mesh.getNumNodes(); n++) {
            gmsh::model::mesh::getNode(nodeTags[n], coord, parametricCoord);
            u_prev[n] = f(coord);
        }

        // Convection vector
        std::vector<double> a = {3, 0, 0};

        // The mass matrix is invariant along iteration
        mesh.precomputeMassMatrix();

        // Element variables
        double elFlux[mesh.getElNumNodes()];
        double elStiffvector[mesh.getElNumNodes()];

        int step = 0;
        for(double t=config.timeStart, tDisplay =0; t<config.timeEnd; t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            u = u_prev;

            // Savings
            if(step==0 || tDisplay>=0.05){
                tDisplay = 0;
                for (int i = 0; i < nodeTags.size(); ++i)
                    solution[i][0] = u[i];
                gmsh::view::addModelData(viewTag, step, names[0], "NodeData", nodeTags, solution, t, 1);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) +
                                    "s] Step number : " + std::to_string(step));
            }

            mesh.precomputeFlux(a.data(), u.data());
            for(int el=0; el<mesh.getElNum(); ++el) {

                mesh.getElFlux(el, elFlux);
                mesh.getElStiffVector(el, a.data(), u.data(), elStiffvector);

                lapack::minus(elStiffvector, elFlux, N);
                lapack::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u_prev[el*N], config.timeStep, N);
            }

        }

        gmsh::view::write(viewTag, "data.msh");
    }
}
