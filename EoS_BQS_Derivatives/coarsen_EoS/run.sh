#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 72:00:00

#./coarsen_EoS ../Thermodynamics_dense/EoS_Taylor_AllMu.dat

# Run with: sbatch --export=ALL,filename=...,r0=... run.sh

./coarsen_EoS $filename $r0

