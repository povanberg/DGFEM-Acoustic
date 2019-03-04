#include "solver.h"
#include <configParser.h>
#include <vector>
#include <Mesh.h>
#include "logger.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <gmsh.h>
#include <Eigen/SparseCholesky>

double integrateGauss(){
    return 0.;
}

namespace solver {

    // Initial confiditon
    double f(Eigen::Vector3d x){
        if(x(0)>1 && x(0)<1.5 && x(1)>0.3 && x(1)<0.8)
            return 1;
        else
            return 0;
    }

    void solveForwardEuler(Mesh &mesh, const Config config) {

        // Solution
        std::vector<int> nodeTags;
        mesh.getNodeTags(nodeTags);
        Eigen::VectorXd u = Eigen::ArrayXd::Zero(nodeTags.size());
        Eigen::VectorXd u_next = Eigen::ArrayXd::Zero(nodeTags.size());

        // Set Initial conditions
        std::vector<double> coord;
        Eigen::VectorXd coordEigen;
        std::vector<double> parametricCoord;
        for(int i=0; i<nodeTags.size(); ++i){
            gmsh::model::mesh::getNode(nodeTags[i], coord, parametricCoord);
            coordEigen = Eigen::Map<Eigen::VectorXd>(coord.data(), coord.size());
            u(i) = f(coordEigen);
        }
        mesh.setData(u);

        // Compute mass matrix
        Eigen::SparseMatrix<double> M;
        mesh.getMassMatrix(M);
        // Inverse
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(M);
        Eigen::SparseMatrix<double> I(M.rows(), M.cols());
        I.setIdentity();
        Eigen::SparseMatrix<double> M_inv = solver.solve(I);

        // Compute stiffness (or convection) matrix
        Eigen::Vector3d a = {1, 1, 0}; // advection coefficient
        Eigen::SparseMatrix<double> K;
        mesh.getStiffMatrix(K, a);

        // Save results example
        int viewTag = gmsh::view::add("TODO_NameConfig");
        std::vector<std::string> names;
        gmsh::model::list(names);
        std::vector<std::vector<double>> solution(nodeTags.size(), std::vector<double>(1));

        // Forward Euler
        double t0 = 0;
        double tf = 1;
        int step = 0;
        double dt = 0.1;
        for(double t=t0; t<tf; t+=dt, ++step) {

            // Savings
            for(int i=0; i<nodeTags.size(); ++i)
                solution[i][0] = u(i);
            gmsh::view::addModelData(viewTag, step, names[0], "NodeData", nodeTags, solution, t, 1);

            // Get Flux
            Eigen::VectorXd F;
            mesh.getFlux(F, a);

            // FE step
            u_next = u + dt*M_inv*(K*u - F);

            // Update previous value: u <- u_next
            u = u_next;
            mesh.setData(u);

            Log("[%f/%f] Simulation time step.", t, tf);
        }

        gmsh::view::write(viewTag, "data.msh");
    }
}