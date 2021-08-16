#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 08:00:00

./coarsen_EoS ../Thermodynamics_dense/EoS_Taylor_AllMu.dat

