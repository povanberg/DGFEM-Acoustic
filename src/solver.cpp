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
        return (1) * exp(-((x[0] - 10) * (x[0] - 10) + (x[1]- 0) * (x[1]- 0)));
    }

    void solveForwardEuler(Mesh &mesh, const Config config) {

        // Solution vectors
        std::vector<int> nodeTags(&mesh.elNodeTag(0), &mesh.elNodeTag(0) + mesh.getNumNodes());
        std::vector<double> u(mesh.getNumNodes());
        std::vector<double> u_prev(mesh.getNumNodes());

        // The mass matrix is invariant along iteration
        mesh.precomputeMassMatrix();

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
        std::vector<double> a = {1, 0, 0};

        // Forward Euler
        double t0 = 0;
        double tf = 0.1;
        int step = 0;
        double dt = 0.001;
        for(double t=t0; t<tf; t+=dt, ++step) {

            gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(tf) +
                                "s] Step number : " + std::to_string(step));

            u = u_prev;

            // Savings
            for (int i = 0; i < nodeTags.size(); ++i)
                solution[i][0] = u[i];
            gmsh::view::addModelData(viewTag, step, names[0], "NodeData", nodeTags, solution, t, 1);

            // Recompute flux for all surfaces
            mesh.precomputeFlux(a.data(), u.data());

            // Compilator is clever and will instantiate out of loop
            double elFlux[mesh.getElNumNodes()];
            double elStiffvector[mesh.getElNumNodes()];

            for(int el=0; el<mesh.getElNum(); ++el) {
                mesh.getElFlux(el, elFlux);
                mesh.getElStiffVector(el, a.data(), u.data(), elStiffvector);

                // u_next = u + dt*M^-1*(K-F)
                /*char TRANS = 'T';
                int INC = 1;
                int N = mesh.getElNumNodes();
                double alpha=dt;
                double beta=1;
                lapack::minus(elStiffvector, elFlux, N);
                dgemv_(TRANS, N, N, alpha, &mesh.elMassMatrix(el), N, &elStiffvector[0], INC, beta, &u_prev[el*N], INC);*/

                // Eigen version
                int N = mesh.getElNumNodes();
                Eigen::Map<Eigen::MatrixXd> _elMassMatrix(&mesh.elMassMatrix(el), N, N);
                Eigen::Map<Eigen::VectorXd> _elFlux(elFlux, N);
                Eigen::Map<Eigen::VectorXd> _elStiffVector(elStiffvector, N);
                Eigen::Map<Eigen::VectorXd> _u(&u_prev[el*N], N);

                _u = _u + dt*_elMassMatrix*(_elStiffVector-_elFlux);

            }

        }

        gmsh::view::write(viewTag, "data.msh");
    }
}
