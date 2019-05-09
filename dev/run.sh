#!/bin/bash
# Submission script for cluster
#SBATCH --job-name=AcousticDGFEM
#SBATCH --time=01:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1024 # megabytes 
#SBATCH --partition=defq 
#SBATCH --output=out.txt

module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=16
export OMP_CANCELLATION=true

cd $HOME/MATH0471-DG/build/bin
srun ./dgalerkin ../../2d/disk.msh ../../config/config.conf
