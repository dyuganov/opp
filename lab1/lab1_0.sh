#!/bin/bash

#PBS -l willtime=00:05:00
#PBS -l select=1:ncpus=1:mem=4000m,place=scatter

cd $PBS_O_WORKDIR

SRC=lab1.cpp
EXE=$SRC.out

echo "Run on node: `uname -n`"
echo "Working dir: $PBS_O_WOKRDIR"

echo "Compiling $SRC"
g++ $SRC -o $EXE

echo "Running $EXE"
./$EXE

