#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>
#include <gmsh.h>
#include <cstdio>
#include <sstream>
#include <utils.h>

#include "configParser.h"

/**
 * Parse et load config file.
 * Update accordingly the global variables.
 */
namespace config{
    std::vector<std::string> split(std::string str, char delimiter) {
        std::vector<std::string> internal;
        std::stringstream ss(str); // Turn the string into a stream.
        std::string tok;

        while(getline(ss, tok, delimiter)) {
            internal.push_back(tok);
        }

        return internal;
    }

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
            config.elementType = configMap["elementType"];
            config.timeIntMethod = configMap["timeIntMethod"];
            config.saveFile = configMap["saveFile"];
            config.numThreads = std::stoi(configMap["numThreads"]);
            config.numThreads = config.numThreads==1? 0 : config.numThreads;
            config.v0[0] = std::stod(configMap["v0_x"]);
            config.v0[1] = std::stod(configMap["v0_y"]);
            config.v0[2] = std::stod(configMap["v0_z"]);
            config.rho0 = std::stod(configMap["rho0"]);
            config.c0 = std::stod(configMap["c0"]);

            for(std::map<std::string, std::string>::iterator iter = configMap.begin(); iter != configMap.end(); ++iter) {
                std::string key = iter->first;
                if(key.find("source") == 0) {
                    std::vector<std::string> sep = split(iter->second, ',');

                    double x = std::stod(sep[1]);
                    double y = std::stod(sep[2]);
                    double z = std::stod(sep[3]);
                    double size = std::stod(sep[4]);
                    double amp = std::stod(sep[5]);
                    double freq = std::stod(sep[6]);
                    double phase = std::stod(sep[7]);
                    double duration = std::stod(sep[8]);
                    double pole;
                    if (sep[0] == "dipole") {
                        pole = 1;
                        std::vector<double> source1 = {pole, x-size, y, z, size/2., amp, freq, phase, duration};
                        std::vector<double> source2 = {pole, x+size, y, z, size/2., amp, freq, phase+M_PI, duration};
                        config.sources.push_back(source1);
                        config.sources.push_back(source2);
                    }
                    else if (sep[0] == "quadrupole") {
                        pole = 2;
                        std::vector<double> source1 = {pole, x-size, y, z, size/2., amp, freq, phase, duration};
                        std::vector<double> source2 = {pole, x+size, y, z, size/2., amp, freq, phase, duration};
                        std::vector<double> source3 = {pole, x, y-size, z, size/2., amp, freq, phase+M_PI, duration};
                        std::vector<double> source4 = {pole, x, y+size, z, size/2., amp, freq, phase+M_PI, duration};
                        config.sources.push_back(source1);
                        config.sources.push_back(source2);
                        config.sources.push_back(source3);
                        config.sources.push_back(source4);
                    }
                    else {
                        pole = 0;
                        std::vector<double> source = {pole, x, y, z, size, amp, freq, phase, duration};
                        config.sources.push_back(source);
                    }
                }
		else if(key.find("initialCondtition") == 0) {
		    std::vector<std::string> sep = split(iter->second, ',');
                    double x = std::stod(sep[1]);
                    double y = std::stod(sep[2]);
                    double z = std::stod(sep[3]);
                    double size = std::stod(sep[4]);
                    double amp = std::stod(sep[5]);
		    std::vector<double> init1 = {0, x, y, z, size, amp};
		    config.initConditions.push_back(init1);
		}
            }

            std::remove(&config.saveFile[0]);

            std::string physName;
            gmsh::vectorpair m_physicalDimTags;
            int bcDim = gmsh::model::getDimension() - 1;
            gmsh::model::getPhysicalGroups(m_physicalDimTags, bcDim);
            for(int p=0; p<m_physicalDimTags.size(); ++p) {
                gmsh::model::getPhysicalName(m_physicalDimTags[p].first, m_physicalDimTags[p].second, physName);
                if (configMap[physName].find("Absorbing") == 0) {
                    config.physBCs[m_physicalDimTags[p].second] = std::make_pair("Absorbing", 0);
                }
                else if (configMap[physName].find("Reflecting") == 0) {
                    config.physBCs[m_physicalDimTags[p].second] = std::make_pair("Reflecting", 0);
                }
                else {
                    gmsh::logger::write("Not specified or supported boundary conditions.");
                }
            }

            gmsh::logger::write("==================================================");
            gmsh::logger::write("Simulation parameters : ");
            gmsh::logger::write("Time step : " + std::to_string(config.timeStep));
            gmsh::logger::write("Final time : " + std::to_string(config.timeEnd));
            gmsh::logger::write("Mean flow velocity : (" + std::to_string(config.v0[0]) + ","
                                + std::to_string(config.v0[1]) + ","
                                + std::to_string(config.v0[2]) + ")");
            gmsh::logger::write("Mean density : " + std::to_string(config.rho0));
            gmsh::logger::write("Speed of sound : " + std::to_string(config.c0));
            gmsh::logger::write("Solver : " + config.timeIntMethod);
        }
        else {
            throw;
        }

        return config;
    }
}


