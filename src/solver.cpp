#include <vector>
#include <gmsh.h>
#include <utils.h>
#include <iostream>
#include "configParser.h"
#include "Mesh.h"

namespace solver {

    // Common variables to all solver
    int elNumNodes;
    std::vector<int> elNodeTags;
    std::vector<double> elFlux;
    std::vector<double> elStiffvector;
    // Gmsh
    int g_viewTag;
    std::vector<std::string> g_names;

    void numStep(Mesh &mesh, Config config, std::vector<double> &u, std::vector<double> &a, double beta){

        mesh.precomputeFlux(a.data(), u.data());
        for(int el=0; el<mesh.getElNum(); ++el){

            mesh.getElFlux(el, elFlux.data());
            mesh.getElStiffVector(el, a.data(), u.data(), elStiffvector.data());

            lapack::minus(elStiffvector.data(), elFlux.data(), elNumNodes);
            lapack::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[el*elNumNodes],
                          config.timeStep, beta, elNumNodes);
        }
    }

    // Solve using forward explicit scheme. O(h)
    // u : initial solution vector    |   mesh   : ...
    // a : convection vector          |   config : ...
    void forwardEuler(std::vector<double> &u, std::vector<double> &a, Mesh &mesh, Config config) {

        std::vector<std::vector<double>> g_u(mesh.getElNodeTags().size(), std::vector<double>(1));
        g_viewTag = gmsh::view::add("Results");
        gmsh::model::list(g_names);

        elNumNodes = mesh.getElNumNodes();
        elNodeTags = mesh.getElNodeTags();
        elFlux.resize(mesh.getElNumNodes());
        elStiffvector.resize(mesh.getElNumNodes());

        mesh.precomputeMassMatrix();
        mesh.setNumFlux(config.flux, a.data(), config.fluxCoeff);

        for(double t=config.timeStart, tDisplay=0, step=0;
                   t<=config.timeEnd;
                   t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            // Savings. (Done at start of step to catch initial configuration)
            if(tDisplay>=config.timeRate || step==0){
                tDisplay = 0;
                for (int n = 0; n < elNodeTags.size(); ++n) {
                    g_u[n][0] = u[n];
                }
                gmsh::view::addModelData(g_viewTag, step, g_names[0], "NodeData", elNodeTags, g_u, t, 1);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) +
                                    "s] Step number : " + std::to_string(step));
            }

            numStep(mesh, config, u, a, 1.0);

            mesh.enforceDiricheletBCs(u.data());
        }

        gmsh::view::write(g_viewTag, "data.msh");
    }

    // Solve using explicit Runge-Kutta integration method. O(h^4)
    // u : initial solution vector    |   mesh   : ...
    // a : convection vector          |   config : ...
    void rungeKutta(std::vector<double> &u, std::vector<double> &a, Mesh &mesh, Config config) {

        std::vector<std::vector<double>> g_u(mesh.getElNodeTags().size(), std::vector<double>(1));
        g_viewTag = gmsh::view::add("Results");
        gmsh::model::list(g_names);

        elNumNodes = mesh.getElNumNodes();
        elNodeTags = mesh.getElNodeTags();
        elFlux.resize(mesh.getElNumNodes());
        elStiffvector.resize(mesh.getElNumNodes());

        std::vector<double> k1;
        std::vector<double> k2;
        std::vector<double> k3;
        std::vector<double> k4;

        mesh.precomputeMassMatrix();
        mesh.setNumFlux(config.flux, a.data(), config.fluxCoeff);

        for(double t=config.timeStart, tDisplay=0, step=0;
            t<=config.timeEnd;
            t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            // Savings. (Done at start of step to catch initial configuration)
            if(tDisplay>=config.timeRate || step==0){
                tDisplay = 0;
                for (int n = 0; n < elNodeTags.size(); ++n) {
                    g_u[n][0] = u[n];
                }
                gmsh::view::addModelData(g_viewTag, step, g_names[0], "NodeData", elNodeTags, g_u, t, 1);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) +
                                    "s] Step number : " + std::to_string(step));
            }

            k1 = k2 = k3 = k4 = u;
            numStep(mesh, config, k1, a, 0); // k1
            lapack::plusTimes(k2.data(), k1.data(), 0.5, k2.size()); // k2
            numStep(mesh, config, k2, a, 0);
            lapack::plusTimes(k3.data(), k2.data(), 0.5, k3.size()); // k3
            numStep(mesh, config, k3, a, 0);
            lapack::plusTimes(k4.data(), k3.data(), 1.0, k4.size()); // k4
            numStep(mesh, config, k4, a, 0);

            for(int i=0; i<u.size(); ++i){
                u[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
            }

            mesh.enforceDiricheletBCs(u.data());
        }

        gmsh::view::write(g_viewTag, "data.msh");
    }

}
