#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>
#include "configParser.h"
#include "logger.h"

// Parse et load config file.
// Update accordingly the global variable.
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

        config.dim = std::stoi(configMap["dim"]);
        log("%s successfully loaded.", name.c_str());
    }
    else {
        error("Couldn't open config file for reading.");
        log("Default values will be used.");
    }

    return config;
}

