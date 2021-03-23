#!/bin/bash

#PBS -l walltime=00:15:00
#PBS -l select=1:ncpus=2:mpiprocs=2:mem=4000m,place=scatter

cd $PBS_O_WORKDIR

MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin
SRC=MPI_v2.cpp
SRCLIB=Matrix.h
EXE=$SRC.out
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')

echo "Run on node: `uname -n`"
echo "Working dir: $PBS_O_WORKDIR"
echo "Proc num: $MPI_NP"

echo "Compiling $SRC $SRCLIB"
$MPI_PATH/mpicxx -std=c++11 $SRC $SRCLIB -o $EXE

echo "Running $EXE"
$MPI_PATH/mpirun -machinefile $PBS_NODEFILE -np $MPI_NP ./$EXE
