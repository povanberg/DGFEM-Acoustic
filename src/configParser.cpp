#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>

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
            config.flux = configMap["flux"];
            config.elementType = configMap["elementType"];
            config.timeIntMethod = configMap["timeIntMethod"];
            config.saveFile = configMap["saveFile"];
        }
        else {
            throw;
        }

        return config;
    }
}


