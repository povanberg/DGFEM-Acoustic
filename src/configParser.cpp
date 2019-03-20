#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>
#include <gmsh.h>

#include "configParser.h"

// Parse et load config file.
// Update accordingly the global variables.
namespace config{
    Config parseConfig(std::string name){

        Config config;

        std::ifstream cFile(name);
        if (cFile.is_open())
        {
            std::string line;
            std::map<std::string, std::string> configMap;
            while(getline(cFile, line)){
                line.erase(std::remove_if(line.begin(), line.end(), isspace),
                           line.end());
                if(line[0] == '#' || line.empty())
                    continue;
                auto delimiterPos = line.find("=");
                auto name = line.substr(0, delimiterPos);
                auto value = line.substr(delimiterPos + 1);
                configMap[name] = value;
            }

            config.timeStart = std::stod(configMap["timeStart"]);
            config.timeEnd = std::stod(configMap["timeEnd"]);
            config.timeStep = std::stod(configMap["timeStep"]);
            config.timeRate = std::stod(configMap["timeRate"]);
            config.flux = configMap["flux"];
            config.fluxCoeff = std::stod(configMap["fluxCoeff"]);
            config.elementType = configMap["elementType"];
            config.timeIntMethod = configMap["timeIntMethod"];
            config.saveFile = configMap["saveFile"];

            // Config physical group must match gmsh physical group
            std::string physName;
            gmsh::vectorpair m_physicalDimTags;
            gmsh::model::getPhysicalGroups(m_physicalDimTags);
            for(int p=0; p<m_physicalDimTags.size(); ++p) {
                gmsh::model::getPhysicalName(m_physicalDimTags[p].first, m_physicalDimTags[p].second, physName);
                if(configMap[physName].find("Dirichelet") == 0) {
                    int valueBC = std::stod(configMap[physName].substr(10));
                    config.physBCs[m_physicalDimTags[p].second] = std::make_pair("Dirichelet", valueBC);
                }
                else if (configMap[physName].find("Neumann") == 0) {
                    int valueBC = std::stod(configMap[physName].substr(7));
                    config.physBCs[m_physicalDimTags[p].second] = std::make_pair("Neumann", valueBC);
                }
                else {
                    gmsh::logger::write("Not supported boundary conditions. Dirichelet 0 applied.");
                    config.physBCs[m_physicalDimTags[p].second] = std::make_pair("Dirichelet", 0);
                }

            }
        }
        else {
            throw;
        }

        return config;
    }
}


