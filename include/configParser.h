#include <string>

#ifndef DGALERKIN_CONFIG_H
#define DGALERKIN_CONFIG_H

struct Config {
    int dim;
    std::string numEquation;
    std::string timeIntMethod;
};

Config parseConfig(std::string name);

#endif
