#ifndef DGALERKIN_LOGGER_H
#define DGALERKIN_LOGGER_H

#define Log(msg, ...) fprintf(stdout,"[INFO] " msg "\n", ##__VA_ARGS__)
#define Error(msg, ...) fprintf(stdout,"[ERR] " msg "\n", ##__VA_ARGS__)

#endif //DGALERKIN_LOGGER_H
