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
    /**
     * The DGarlekin solver requires 2 arguments
     * 1 : the Mesh file (.msh)
     * 2 : the config file (.conf)
     *
     * e.g. ./dgarlerkin mymesh.msh myconfig.conf
     */
    if(argc!=3){ return E2BIG; }
    std::string msh_name = argv[1];
    std::string config_name = argv[2];

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(msh_name);

    Config config = config::parseConfig(config_name);
    gmsh::logger::write("Config loaded : " + config_name);

    Mesh mesh(msh_name, config);

    /**
     * Initialize the solution:
     * Here you can impose the initial condition
     */
    std::vector<std::vector<double>> u(4,std::vector<double>(mesh.getNumNodes()));
    for(int n=0; n<mesh.getNumNodes(); n++){
        std::vector<double> coord, paramCoord;
        gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord);
        u[0][n] = exp(-((coord[0] - 0) * (coord[0] - 0) +
                        (coord[1]- 0) * (coord[1] - 0) +
                        (coord[2]- 0) * (coord[2]- 0))/1);
        //u[0][n] = 0;
        u[1][n] = 0;
        u[2][n] = 0;
        u[3][n] = 0;
    }

    if(config.timeIntMethod == "Euler1")
        solver::forwardEuler(u, mesh, config);
    else if(config.timeIntMethod == "Runge-Kutta")
        solver::rungeKutta(u, mesh, config);

    gmsh::finalize();

    return EXIT_SUCCESS;
}
