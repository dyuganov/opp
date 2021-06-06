#!/bin/bash

#PBS -l walltime=00:05:00
#PBS -l select=2:ncpus=8:mpiprocs=8:mem=4000m,place=scatter

cd $PBS_O_WORKDIR

MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin
SRC=MPI_v2_2.cpp
SRCLIB=Matrix.h
EXE=$SRC.out
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')

echo "Run on node: `uname -n`"
echo "Number of MPI process: $MPI_NP"
echo "Workong dir: $PBS_O_WORKDIR"

echo "Compiling $SRC $SRCLIB"
mpiicc -std=c++11 $SRC $SRLIB -o $EXE
# $MPI_PATH/mpicxx -std=c++11 $SRC $SRLIB -o $EXE

echo "Running $SRC $SRCLIB"
$MPI_PATH/mpirun -trace -machinefile $PBS_NODEFILE -np $MPI_NP ./$EXE
