#include <vector>
#include <gmsh.h>
#include <utils.h>
#include <iostream>
#include "configParser.h"
#include "Mesh.h"

namespace solver {

    // Initial condititon
    double f(std::vector<double> x) {
        //return (1) * exp(-((x[0] + 3.2) * (x[0] + 3.2))/(0.5));
        return (1) * exp(-((x[0] - 10) * (x[0] - 10) + (x[1]- 0) * (x[1]- 0))/0.5);
    }

    void dtfu(Mesh &mesh, Config config, std::vector<double> &u_prev, std::vector<double> &a, double beta, int N){

        // Element variables
        double elFlux[mesh.getElNumNodes()];
        double elStiffvector[mesh.getElNumNodes()];

        mesh.precomputeFlux(a.data(), u_prev.data());
        for(int el=0; el<mesh.getElNum(); ++el){
                    
            mesh.getElFlux(el, elFlux);
            mesh.getElStiffVector(el, a.data(), u_prev.data(), elStiffvector);
            lapack::minus(elStiffvector, elFlux, N);
            lapack::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u_prev[el*N], config.timeStep, beta, N);
        }
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
        std::vector<double> a = {1, 1, 0};

        // The mass matrix is invariant along iteration
        mesh.precomputeMassMatrix();

        int step = 0;
        for(double t=config.timeStart, tDisplay =0; t<config.timeEnd; t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            u = u_prev;

            // Savings
            if(step==0 || tDisplay>=config.timeRate){
                tDisplay = 0;
                for (int i = 0; i < nodeTags.size(); ++i)
                    solution[i][0] = u[i];
                gmsh::view::addModelData(viewTag, step, names[0], "NodeData", nodeTags, solution, t, 1);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) +
                                    "s] Step number : " + std::to_string(step));
            }

            dtfu(mesh, config, u_prev, a, 1.0, N);
        }

        gmsh::view::write(viewTag, "data.msh");
    }

    void solveRungeKutta(Mesh &mesh, Config config) {

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
        std::vector<double> a = {1, 1, 0};

        // The mass matrix is invariant along iteration
        mesh.precomputeMassMatrix();

        int step = 0;
        for(double t=config.timeStart, tDisplay =0; t<config.timeEnd; t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            u = u_prev;

            // Savings
            if(step==0 || tDisplay>=config.timeRate){
                tDisplay = 0;
                for (int i = 0; i < nodeTags.size(); ++i)
                    solution[i][0] = u[i];
                gmsh::view::addModelData(viewTag, step, names[0], "NodeData", nodeTags, solution, t, 1);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) +
                                    "s] Step number : " + std::to_string(step));
            }

            // k vectors
            std::vector<double> k1 = u_prev;
            std::vector<double> k2 = u_prev;
            std::vector<double> k3 = u_prev;
            std::vector<double> k4 = u_prev;

            dtfu(mesh, config, k1, a, 0, N); // k1
            lapack::plusTimes(k2.data(), k1.data(), 0.5, k2.size()); // k2
            dtfu(mesh, config, k2, a, 0, N);
            lapack::plusTimes(k3.data(), k2.data(), 0.5, k3.size()); // k3
            dtfu(mesh, config, k3, a, 0, N);
            lapack::plusTimes(k4.data(), k3.data(), 1.0, k4.size()); // k4
            dtfu(mesh, config, k4, a, 0, N);

            for(int i=0; i<u_prev.size(); ++i){
                u_prev[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
            }
        }

        gmsh::view::write(viewTag, "data.msh");

    }
}
