#!/bin/sh

echo "******************************************";
echo "*     Discontinuous galerkin setup.      *";
echo "******************************************";
echo

echo "[0] Create some project directories if don't exist.";

mkdir 3rdParty

mkdir results

echo "[1] Get depedencies and external libraries.";


dpkg -s cmake > /dev/null 2>&1;
if [ $? -eq 0 ]; then
	echo "mCmake found.";
else
	echo "Cmake not found, installing...";
	sudo apt-get -y install cmake;
	echo "Cmake installed.";
fi

dpkg -s g++ > /dev/null 2>&1;
if [ $? -eq 0 ]; then
	echo "G++ found.";
else
	echo "G++ not found, installing...";
	sudo apt-get -y install g++;
	export CC=gcc;
	export CXX=g++;
	echo "G++ installed.";
fi

dpkg -s gfortran > /dev/null 2>&1;
if [ $? -eq 0 ]; then
	echo "Gfortran found.";
else
	echo "Gfortran not found, installing...";
	sudo apt-get -y install gfortran;
	echo "Gfortran installed.";
fi

dpkg -s libblas-dev liblapack-dev > /dev/null 2>&1;
if [ $? -eq 0 ]; then
	echo "Lapack/Blas found.";
else
	echo "Lapack/Blas not found, installing...";
	sudo apt-get -y install libblas-dev liblapack-dev;
	echo "Lapack/Blas installed.";
fi

#if [ ! -d "3rdParty/gmsh" ]; then
#	echo "Gmsh not found, installing...";
#	wget http://gmsh.info/bin/Linux/gmsh-4.1.5-Linux64-sdk.tgz
#	tar -xf gmsh-4.1.5-Linux64-sdk.tgz
#	rm -rf gmsh-4.1.5-Linux64-sdk.tgz
#	mv gmsh-4.1.5-Linux64-sdk 3rdParty/gmsh
#	echo "Gmsh installed."
#else
#       echo "Gmsh found.";

if [ ! -d "3rdParty/gmsh" ]; then
	echo "Gmsh not found, installing...";
	wget http://gmsh.info/bin/Linux/gmsh-4.10.5-Linux64-sdk.tgz
	tar -xf gmsh-4.10.5-Linux64-sdk.tgz
	rm -rf gmsh-4.10.5-Linux64-sdk.tgz
	mv gmsh-4.10.5-Linux64-sdk 3rdParty/gmsh
	echo "Gmsh installed."
else
       echo "Gmsh found.";
       
fi
cd 3rdParty/gmsh/
export FC=gfortran
export PATH=${PWD}/bin:${PWD}/lib:${PATH}
export INCLUDE=${PWD}/include:${INCLUDE}
export LIB=${PWD}/lib:${LIB}
export PYTHONPATH=${PWD}/lib:${PYTHONPATH} 
export DYLD_LIBRARY_PATH=${PWD}/lib:${DYLD_LIBRARY_PATH}
cd ../../


if [ ! -d "3rdParty/eigen" ]; then
	echo "Eigen not found, installing...";
  	wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
 	tar -xf eigen-3.4.0.tar.gz
 	rm -rf eigen-3.4.0.tar.gz
 	mv eigen-3.4.0 3rdParty/eigen
	echo "Eigen installed."
else
	echo "Eigen found."
fi
cd 3rdParty/eigen/
export INCLUDE=${PWD}:${INCLUDE}
cd ../../


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
