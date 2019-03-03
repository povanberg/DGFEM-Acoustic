# This file setup an environment which allows you to compile the code
# on Linux or macOS using the default compiler (gcc or clang)
#
# HOW TO USE THIS FILE?
#   open a terminal
#   . ./devenv-linux-macos.sh     # the first dot is important!
#   mkdir build
#   cd build
#   cmake ..
#   make
#   [executables are built in the bin/ folder]
#   ctest

# set the location of gmsh SDK ( **MODIFY THIS LINE FOR YOUR SYSTEM** )
GMSHSDK=~/Desktop/gmsh-sdk
# location of eigen
EIGEN=~/Desktop/eigen

# where are gmsh and gmsh-**.so ?
export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:${PATH}
# where is gmsh.h ?
export INCLUDE=${GMSHSDK}/include:${INCLUDE}
# eigen library
export INCLUDE=${EIGEN}:${INCLUDE}
# where is gmsh.lib ?
export LIB=${GMSHSDK}/lib:${LIB}
# where is gmsh.py ? (required only if you want to use the python API)
export PYTHONPATH=${GMSHSDK}/lib:${PYTHONPATH}
# the following command is only useful for macOS 
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:${DYLD_LIBRARY_PATH}
