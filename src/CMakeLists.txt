# programs

SET(SRCS dgalerkin.cpp configParser.cpp ../include/configParser.h Mesh.cpp ../include/Mesh.h utils.cpp ../include/utils.h solver.cpp ../include/solver.h)
ADD_EXECUTABLE(dgalerkin ${SRCS})
TARGET_LINK_LIBRARIES(dgalerkin ${GMSH_LIBRARIES} ${LAPACKBLAS_LIBRARIES})

SET(SRCS testUnit.cpp configParser.cpp ../include/configParser.h Mesh.cpp ../include/Mesh.h utils.cpp ../include/utils.h solver.cpp ../include/solver.h)
ADD_EXECUTABLE(testUnit ${SRCS})
TARGET_LINK_LIBRARIES(testUnit ${GMSH_LIBRARIES} ${LAPACKBLAS_LIBRARIES})
