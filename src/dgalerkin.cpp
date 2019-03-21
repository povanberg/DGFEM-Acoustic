#include <cstdio>
#include <gmsh.h>
#include <errno.h>
#include <iostream>

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

    // Retrieve parameters required to run the
    // Discontinuous Galerkin simulation.
    Mesh mesh(msh_name, config);

    // Convection vector
    std::vector<double> a = {3, 3, 0};

    // Initialize the solution
    std::vector<double> u(mesh.getNumNodes());
    for(int n=0; n<mesh.getNumNodes(); n++) {
        std::vector<double> coord, paramCoord;
        gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord);
        // Gaussian
        u[n] = exp(-((coord[0] - 10) * (coord[0] - 10) + (coord[1]- 0) * (coord[1]- 0))/0.5);
    }

    if(config.timeIntMethod == "Euler1")
        solver::forwardEuler(u, a, mesh, config);
    else if(config.timeIntMethod == "Runge-Kutta")
        solver::rungeKutta(u, a, mesh, config);

    gmsh::finalize();

    return EXIT_SUCCESS;
}
