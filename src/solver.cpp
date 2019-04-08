#include <vector>
#include <gmsh.h>
#include <utils.h>
#include <iostream>
#include <chrono>
#include <omp.h>
#include <sstream>

#include "configParser.h"
#include "Mesh.h"

namespace solver {

    // Common variables to all solver
    int elNumNodes;
    std::vector<int> elTags;
    std::vector<double> elFlux;
    std::vector<double> elStiffvector;
    // Gmsh
    int g_viewTag;
    std::vector<std::string> g_names;

    void numStep(Mesh &mesh, Config config, std::vector<double> &u, std::vector<double> &a, double beta){

        mesh.precomputeFlux(a.data(), u.data());

        #pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector) num_threads(config.numThreads)
        for(int el=0; el<mesh.getElNum(); ++el){

            mesh.getElFlux(el, elFlux.data());
            mesh.getElStiffVector(el, a.data(), u.data(), elStiffvector.data());

            eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes);
            eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[el*elNumNodes],
                        config.timeStep, beta, elNumNodes);
        }
    }

    // Solve using forward explicit scheme. O(h)
    // u : initial solution vector    |   mesh   : ...
    // a : convection vector          |   config : ...
    void forwardEuler(std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &a, Mesh &mesh, Config config) {

        elNumNodes = mesh.getElNumNodes();
        elTags = std::vector<int>(&mesh.elTag(0), &mesh.elTag(0)+mesh.getElNum());
        elFlux.resize(mesh.getElNumNodes());
        elStiffvector.resize(mesh.getElNumNodes());

        g_viewTag = gmsh::view::add("Results");
        gmsh::model::list(g_names);
        std::vector<std::vector<double>> g_u(mesh.getElNum(), std::vector<double>(elNumNodes));

        mesh.precomputeMassMatrix();

        auto start = std::chrono::system_clock::now();

        for(double t=config.timeStart, tDisplay=0, step=0;
                   t<=config.timeEnd;
                   t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            // Savings. (Done at start of step to catch initial configuration)
            if(tDisplay>=config.timeRate || step==0){
                tDisplay = 0;
                for(int el = 0; el < mesh.getElNum(); ++el) {
                    for(int n = 0; n < mesh.getElNumNodes(); ++n) {
                        g_u[el][n] = u[0][el*elNumNodes+n]+u[1][el*elNumNodes+n]+u[2][el*elNumNodes+n]+u[3][el*elNumNodes+n];
                    }
                }
                gmsh::view::addModelData(g_viewTag, step, g_names[0], "ElementNodeData", elTags, g_u, t, 1);
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) + "s] Step number : "
                                    + std::to_string((int) step) + ", Elapsed time: " + std::to_string(elapsed.count()) +"s");
            }

            for(int v=0; v<u.size(); ++v){
                mesh.setNumFlux(config.flux, a[v], config.fluxCoeff);
                numStep(mesh, config, u[v], a[v], 1.0);
                mesh.enforceDiricheletBCs(u[v]);
            }
        }

        gmsh::view::write(g_viewTag, "data.msh");
    }

    // Solve using explicit Runge-Kutta integration method. O(h^4)
    // u : initial solution vector    |   mesh   : ...
    // a : convection vector          |   config : ...
    void rungeKutta(std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &a, Mesh &mesh, Config config) {

        elNumNodes = mesh.getElNumNodes();
        elTags = std::vector<int>(&mesh.elTag(0), &mesh.elTag(0)+mesh.getElNum());
        elFlux.resize(mesh.getElNumNodes());
        elStiffvector.resize(mesh.getElNumNodes());

        g_viewTag = gmsh::view::add("Results");
        gmsh::model::list(g_names);
        std::vector<std::vector<double>> g_u(mesh.getElNum(), std::vector<double>(elNumNodes));

        std::vector<double> k1;
        std::vector<double> k2;
        std::vector<double> k3;
        std::vector<double> k4;

        mesh.precomputeMassMatrix();

        auto start = std::chrono::system_clock::now();

        for(double t=config.timeStart, tDisplay=0, step=0;
            t<=config.timeEnd;
            t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            // Savings. (Done at start of step to catch initial configuration)
            if(tDisplay>=config.timeRate || step==0){
                tDisplay = 0;
                for(int el = 0; el < mesh.getElNum(); ++el) {
                    for(int n = 0; n < mesh.getElNumNodes(); ++n) {
                        g_u[el][n] = u[0][el*elNumNodes+n]+u[1][el*elNumNodes+n]+u[2][el*elNumNodes+n]+u[3][el*elNumNodes+n];
                    }
                }
                gmsh::view::addModelData(g_viewTag, step, g_names[0], "ElementNodeData", elTags, g_u, t, 1);
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) + "s] Step number : "
                                    + std::to_string((int) step) + ", Elapsed time: " + std::to_string(elapsed.count()) +"s");
            }

            for(int v=0; v<u.size(); ++v){
                mesh.setNumFlux(config.flux, a[v], config.fluxCoeff);

                k1 = k2 = k3 = k4 = u[v];
                numStep(mesh, config, k1, a[v], 0); // k1
                eigen::plusTimes(k2.data(), k1.data(), 0.5, k2.size());
                numStep(mesh, config, k2, a[v], 0); // k2
                eigen::plusTimes(k3.data(), k2.data(), 0.5, k3.size());
                numStep(mesh, config, k3, a[v], 0); // k3
                eigen::plusTimes(k4.data(), k3.data(), 1.0, k4.size());
                numStep(mesh, config, k4, a[v], 0); // k4

                for(int i=0; i<u[v].size(); ++i){
                    u[v][i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
                }
                mesh.enforceDiricheletBCs(u[v]);
            }
        }

        gmsh::view::write(g_viewTag, "data.msh");
    }

}
