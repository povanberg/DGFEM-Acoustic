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
    int numNodes;
    std::vector<int> elTags;
    std::vector<double> elFlux;
    std::vector<double> elStiffvector;
    std::vector<std::vector<std::vector<double>>> Flux;
    // Gmsh
    int g_viewTag;
    std::vector<std::string> g_names;

    void numStep(Mesh &mesh, Config config, std::vector<double> &u, std::vector<std::vector<double>> &Flux, double beta, int eq){

        mesh.precomputeFlux(u, Flux, eq);

        #pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector) num_threads(config.numThreads)
        for(int el=0; el<mesh.getElNum(); ++el){

            mesh.getElFlux(el, elFlux.data());
            mesh.getElStiffVector(el, Flux, u, elStiffvector.data());
            eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes);
            eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[el*elNumNodes],
                        config.timeStep, beta, elNumNodes);
        }
    }

    // Solve using forward explicit scheme. O(h)
    // u : initial solution vector    |   mesh   : ...
    // a : convection vector          |   config : ...
    void forwardEuler(std::vector<std::vector<double>> &u, std::vector<std::vector<std::vector<double>>> &a, Mesh &mesh, Config config) {

        /*elNumNodes = mesh.getElNumNodes();
        Numnodes = mesh.getNumNodes();
        elTags = std::vector<int>(&mesh.elTag(0), &mesh.elTag(0)+mesh.getElNum());
        elFlux.resize(elNumNodes);
        elStiffvector.resize(elNumNodes);

        g_viewTag = gmsh::view::add("Results");
        gmsh::model::list(g_names);
        std::vector<std::vector<double>> g_u(mesh.getElNum(), std::vector<double>(elNumNodes));

        mesh.precomputeMassMatrix();
        mesh.setNumFlux(config.flux, config.fluxCoeff);

        auto start = std::chrono::system_clock::now();

        for(double t=config.timeStart, tDisplay=0, step=0;
                   t<=config.timeEnd;
                   t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            // Savings. (Done at start of step to catch initial configuration)
            if(tDisplay>=config.timeRate || step==0){
                tDisplay = 0;
                for(int el = 0; el < mesh.getElNum(); ++el) {
                    for(int n = 0; n < mesh.getElNumNodes(); ++n) {
                        g_u[el][n] = u[0][el*elNumNodes+n];
                    }
                }
                gmsh::view::addModelData(g_viewTag, step, g_names[0], "ElementNodeData", elTags, g_u, t, 1);
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) + "s] Step number : "
                                    + std::to_string((int) step) + ", Elapsed time: " + std::to_string(elapsed.count()) +"s");
            }
            for(int v=0; v<u.size(); ++v){
                numStep(mesh, config, u[v], a[v], 1.0);
            }
            mesh.enforceDiricheletBCs(u);
            mesh.updateFlux(a, u);
        }
        gmsh::view::write(g_viewTag, "data.msh");*/
    }

    // Solve using explicit Runge-Kutta integration method. O(h^4)
    // u : initial solution vector    |   mesh   : ...
    // a : convection vector          |   config : ...
    void rungeKutta(std::vector<std::vector<double>> &u, Mesh &mesh, Config config) {

        // Memory allocation
        elNumNodes = mesh.getElNumNodes();
        numNodes = mesh.getNumNodes();
        elTags = std::vector<int>(&mesh.elTag(0), &mesh.elTag(0)+mesh.getElNum());
        elFlux.resize(elNumNodes);
        elStiffvector.resize(elNumNodes);
        std::vector<std::vector<double>> k1, k2, k3, k4;
        Flux = std::vector<std::vector<std::vector<double>>>(4,
               std::vector<std::vector<double>>(mesh.getNumNodes(),
               std::vector<double>(3)));

        // Gmsh save
        g_viewTag = gmsh::view::add("Results");
        gmsh::model::list(g_names);
        std::vector<std::vector<double>> g_u(mesh.getElNum(), std::vector<double>(elNumNodes));

        // Precomputation
        mesh.precomputeMassMatrix();

        // Time iteration
        auto start = std::chrono::system_clock::now();
        for(double t=config.timeStart, tDisplay=0, step=0;
            t<=config.timeEnd;
            t+=config.timeStep, tDisplay+=config.timeStep, ++step) {

            // Savings. (Done at start of step to catch initial configuration)
            if(tDisplay>=config.timeRate || step==0){
                tDisplay = 0;
                for(int el = 0; el < mesh.getElNum(); ++el) {
                    for(int n = 0; n < mesh.getElNumNodes(); ++n) {
                        g_u[el][n] = u[0][el*elNumNodes+n];
                    }
                }
                gmsh::view::addModelData(g_viewTag, step, g_names[0], "ElementNodeData", elTags, g_u, t, 1);
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) + "s] Step number : "
                                    + std::to_string((int) step) + ", Elapsed time: " + std::to_string(elapsed.count()) +"s");
            }

            // Source
            if(t<5) {
                double s_freq = 4;
                std::vector<double> s_coord = {5, 5, 0};
                std::vector<double> s_size = {0.2, 0.2, 0.1};
                for(int n=0; n<mesh.getNumNodes(); n++) {
                    std::vector<double> coord, paramCoord;
                    gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord);
                    if(abs(coord[0]-s_coord[0])<s_size[0] &&
                       abs(coord[1]-s_coord[1])<s_size[1]) {
                        u[0][n] = sin(6.28*s_freq*t);
                    }
                }
            }

            // Fourth order Runge-Kutta algorithm
            k1 = k2 = k3 = k4 = u;
            mesh.updateFlux(k1, Flux, config.v0, config.c0, config.rho0);
            for(int v=0; v<u.size(); ++v){
                numStep(mesh, config, k1[v], Flux[v], 0, v); // k1
                eigen::plusTimes(k2[v].data(), k1[v].data(), 0.5, numNodes);
            }
            mesh.updateFlux(k2, Flux, config.v0, config.c0, config.rho0);
            for(int v=0; v<u.size(); ++v){
                numStep(mesh, config, k2[v], Flux[v], 0, v); // k2
                eigen::plusTimes(k3[v].data(), k2[v].data(), 0.5, numNodes);
            }
            mesh.updateFlux(k3, Flux, config.v0, config.c0, config.rho0);
            for(int v=0; v<u.size(); ++v){
                numStep(mesh, config, k3[v], Flux[v], 0, v); // k3
                eigen::plusTimes(k4[v].data(), k3[v].data(), 1.0, numNodes);
            }
            mesh.updateFlux(k4, Flux, config.v0, config.c0, config.rho0);
            for(int v=0; v<u.size(); ++v){
                numStep(mesh, config, k4[v], Flux[v], 0, v); // k4
                for(int i=0; i<numNodes; ++i){
                    u[v][i] += (k1[v][i]+2*k2[v][i]+2*k3[v][i]+k4[v][i])/6.0;
                }
            }

        }
        gmsh::view::write(g_viewTag, config.saveFile);
    }
}
