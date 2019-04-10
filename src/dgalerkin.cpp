#include <cstdio>
#include <gmsh.h>
#include <errno.h>
#include <iostream>
#include <omp.h>

#include "Mesh.h"
#include "solver.h"
#include "configParser.h"

int main(int argc, char **argv)
{
    // The DGarlekin solver requires 2 arguments
    // 1: the Mesh file (.msh)
    // 2: the config file (see, configParser.cpp)
    // e.g. ./dgarlerkin mymesh.msh myconfig.conf
    // ------------------------------------------
    if(argc!=3){ return E2BIG; }
    std::string msh_name = argv[1];
    std::string config_name = argv[2];

    // Init Gmsh
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(msh_name);

    // Load config parameters
    Config config = config::parseConfig(config_name);
    gmsh::logger::write("Config loaded : " + config_name);

    // Retrieve parameters required to run the Discontinuous Galerkin simulation.
    Mesh mesh(msh_name, config);

    // Physical flux vector a[variable][node][dimension]
    std::vector<std::vector<std::vector<double>>> a(4,std::vector<std::vector<double>>(mesh.getNumNodes(),std::vector<double>(3)));

    // Initialize the solution
    std::vector<std::vector<double>> u(4,std::vector<double>(mesh.getNumNodes()));
    for(int v=0; v<u.size(); v++){
        for(int n=0; n<mesh.getNumNodes(); n++){
            std::vector<double> coord, paramCoord;
            gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord);
            // Gaussian for P and 0 for v(x,y,z)
            u[0][n] = exp(-((coord[0] - 10) * (coord[0] - 10) + (coord[1]+ 0) * (coord[1]- 0) + (coord[2]- 0) * (coord[2]- 0))/1);
            u[1][n] = 0;
            u[2][n] = 0;
            u[3][n] = 0;
        }
    }

    // Initialise the physical flux vector
    mesh.updateFlux(a, u);

    // Solver
    if(config.timeIntMethod == "Euler1")
        solver::forwardEuler(u, a, mesh, config);
    else if(config.timeIntMethod == "Runge-Kutta")
        solver::rungeKutta(u, a, mesh, config);

    gmsh::finalize();

    return EXIT_SUCCESS;
}
