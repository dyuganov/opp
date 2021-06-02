#!/bin/bash

#PBS -l walltime=00:05:00
#PBS -l select=2:ncpus=8:mpiprocs=8:mem=8000m,place=free

cd $PBS_O_WORKDIR

MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin
SRC=main.cpp
EXE=$SRC.out
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')

echo "Working dir: $PBS_O_WORKDIR"
echo "Proc num: $MPI_NP"

echo "Compiling $SRC"
$MPI_PATH/mpiicc -std=c++11 $SRC -o $EXE

echo "Running $EXE"
$MPI_PATH/mpirun -trace -machinefile $PBS_NODEFILE -np $MPI_NP ./$EXE
