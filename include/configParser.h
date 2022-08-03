#include <fstream>
#include <map>
#include <string>

#ifndef DGALERKIN_CONFIG_H
#define DGALERKIN_CONFIG_H

//! ////////////////////////////////////////////////////////////////
struct Config
{
    // Initial, final time and time step(t>0)
    double timeStart = 0;
    double timeEnd = 1;
    double timeStep = 0.1;
    double timeRate = 0.1;

    // Element Type:
    std::string elementType = "Lagrange";

    // Time integration method
    std::string timeIntMethod = "Euler1";

    // Boundary condition
    // key : physical group Tag
    // value : tuple<BCType, BCValue>
    std::map<int, std::pair<std::string, double>> physBCs;

    // Mean flow parameters
    std::vector<double> v0 = {0, 0, 0};
    double rho0 = 1;
    double c0 = 1;

    // Number of threads
    int numThreads = 1;

    // Sources
    std::vector<std::vector<double>> sources;

    // Initial conditions
    std::vector<std::vector<double>> initConditions;

    // Save file
    std::string saveFile = "results.msh";
};

namespace config
{
    Config parseConfig(std::string name);
}

#endif
