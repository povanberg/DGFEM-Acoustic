# programs

SET(SRCS dgalerkin.cpp configParser.cpp ../include/configParser.h Mesh.cpp ../include/Mesh.h utils.cpp ../include/utils.h solver.cpp ../include/solver.h)
ADD_EXECUTABLE(dgalerkin ${SRCS})
TARGET_LINK_LIBRARIES(dgalerkin ${GMSH_LIBRARIES} ${LAPACKBLAS_LIBRARIES})
