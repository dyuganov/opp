#!/bin/bash

#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=4000m,place=scatter

cd $PBS_O_WORKDIR

MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin
SRC=nonMPI.cpp
EXE=$SRC.out
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')

echo "Run on node: `uname -n`"
echo "Number of MPI process: $MPI_NP"
echo "Workong dir: $PBS_O_WORKDIR"

echo "Compiling $SRC"
$MPI_PATH/mpicxx -std=c++11 $SRC -o $EXE

echo "Running $EXE"
$MPI_PATH/mpirun $PBS_NODEFILE -np $MPI_NP ./$EXE

