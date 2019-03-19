#include <string>

#ifndef DGALERKIN_CONFIG_H
#define DGALERKIN_CONFIG_H

struct Config {
    // Initial, final time and time step(t>0)
    double timeStart = 0;
    double timeEnd =  1;
    double timeStep = 0.1;
    double timeRate = 0.1;

    // Specify the numerical flux
    std::string flux = "average";
    double fluxCoeff=1;

    // Element Type:
    std::string elementType = "Lagrange";

    // Time integration method
    std::string timeIntMethod = "Euler1";

    // Save file
    std::string saveFile = "results.msh";
};

namespace config{
    Config parseConfig(std::string name);
}

#endif
