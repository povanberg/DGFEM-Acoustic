#!/bin/sh

echo "******************************************";
echo "*     Discontinuous galerkin setup.      *";
echo "******************************************";
echo

echo "[1] Get depedencies and external libraries.";

module load cmake/3.11.1
module load gcc/4.9.2

dpkg -s cmake > /dev/null 2>&1;
if [ $? -eq 0 ]; then
	echo "mCmake found.";
else
	echo "Cmake not found, installing...";
	apt-get install cmake;
	echo "Cmake installed.";
fi

dpkg -s g++ > /dev/null 2>&1;
if [ $? -eq 0 ]; then
	echo "G++ found.";
else
	echo "G++ not found, installing...";
	apt-get install g++;
	export CC=gcc;
	export CXX=g++;
	echo "G++ installed.";
fi

dpkg -s gfortran > /dev/null 2>&1;
if [ $? -eq 0 ]; then
	echo "Gfortran found.";
else
	echo "Gfortran not found, installing...";
	apt-get install gfortran;
	echo "Gfortran installed.";
fi

dpkg -s libblas-dev liblapack-dev > /dev/null 2>&1;
if [ $? -eq 0 ]; then
	echo "Lapack/Blas found.";
else
	echo "Lapack/Blas not found, installing...";
	apt-get install libblas-dev liblapack-dev;
	echo "Lapack/Blas installed.";
fi

if [ ! -d "gmsh-4.1.4-Linux64-sdk" ]; then
	echo "Gmsh not found, installing...";
	wget http://gmsh.info/bin/Linux/gmsh-4.1.4-Linux64-sdk.tgz
	tar -xf gmsh-4.1.4-Linux64-sdk.tgz
	rm -rf gmsh-4.1.4-Linux64-sdk.tgz
	echo "Gmsh installed."
else
       echo "Gmsh found.";
fi
cd gmsh-4.1.4-Linux64-sdk/
export FC=gfortran
export PATH=${PWD}/bin:${PWD}/lib:${PATH}
export INCLUDE=${PWD}/include:${INCLUDE}
export LIB=${PWD}/lib:${LIB}
export PYTHONPATH=${PWD}/lib:${PYTHONPATH} 
export DYLD_LIBRARY_PATH=${PWD}/lib:${DYLD_LIBRARY_PATH}
cd ../


if [ ! -d "eigen-eigen-323c052e1731" ]; then
	echo "Eigen not found, installing...";
  	wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
 	tar -xf 3.3.7.tar.gz
 	rm -rf 3.3.7.tar.gz
	echo "Eigen installed."
else
	echo "Eigen found."
fi
cd eigen-eigen-323c052e1731/
export INCLUDE=${PWD}:${INCLUDE}
cd ../


echo "[2] Build sources.";

rm -rf build/  
mkdir build

cd build/
cmake ../ -DCMAKE_BUILD_TYPE=Release  -G "Unix Makefiles" 
make -j4
if [ $? -eq 0 ]; then
    	echo "[end] Everything went successfully.";
else
	echo "[end] Error!";
fi
