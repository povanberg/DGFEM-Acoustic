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

    solver::solveTimeIntegration(mesh, config);

    //gmsh::fltk::run();
    gmsh::finalize();

    return EXIT_SUCCESS;
}
