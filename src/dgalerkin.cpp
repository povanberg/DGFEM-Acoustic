#include <cstdio>
#include <gmsh.h>
#include <errno.h>
#include <iostream>
#include "configParser.h"
#include "logger.h"
#include "Mesh.h"
#include "solver.h"
#include "mass.h"

int main(int argc, char **argv)
{
    // The DGarlekin solver requires 2 arguments
    // 1: the Mesh file (.msh)
    // 2: the config file (see, configParser.cpp)
    // e.g. ./dgarlerkin mymesh.msh myconfig.conf
    if(argc!=3)
    {
        Error("Argument list too long.");
        return E2BIG;
    }

    std::string msh_name = argv[1];
    std::string config_name = argv[2];

    // Load simulation parameters
    /*Config config = parseConfig(config_name);*/

    // Load Mesh
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(msh_name);

    // Create object wrapper (keep in memory variables)
    // Structure: Mesh > Elements > faces > node
    //--------------------------------------------------
    // .msh files must contains physical regions:
    // 'domain' = where the simulation must be run [dim]
    // 'dirichelet' = dirichelet BCs region [dim-1]
    // 'neumann" = neumann BCs region [dim-1] //todo: no yet implemented
    // 'robin" = robin BCs region [dim-1] //todo: no yet implemented
    //--------------------------------------------------
    Mesh mesh(msh_name);

    // Simulation
    /*solver::solveForwardEuler(mesh, config);*/

    // Creates the Mass matrix for each elements in the mesh
    mass::createM(mesh);

    //gmsh::fltk::run();
    gmsh::finalize();

    std::cout << "SUCCESS" << std::endl;

    return EXIT_SUCCESS;
}
