# programs
SET(SRCS dgalerkin.cpp configParser.cpp ../include/configParser.h ../include/logger.h Mesh.cpp ../include/Mesh.h gmshUtils.cpp ../include/gmshUtils.h solver.cpp ../include/solver.h Element.cpp ../include/Element.h Face.cpp ../include/Face.h utils.cpp ../include/utils.h mass.cpp ../include/mass.h .. )
ADD_EXECUTABLE(dgalerkin ${SRCS})
TARGET_LINK_LIBRARIES(dgalerkin ${GMSH_LIBRARIES})