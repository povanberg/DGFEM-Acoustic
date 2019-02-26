#ifndef DGALERKIN_LOGGER_H
#define DGALERKIN_LOGGER_H

#define log(msg, ...) fprintf(stdout,"[INFO] " msg "\n", ##__VA_ARGS__)
#define error(msg, ...) fprintf(stdout,"[ERR] " msg "\n", ##__VA_ARGS__)

#endif //DGALERKIN_LOGGER_H
