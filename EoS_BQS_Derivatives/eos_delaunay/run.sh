#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 72:00:00

./interpolate_ebsq ../Thermodynamics_improved_dense/EoS_Taylor_AllMu.dat ../Thermodynamics_staggered/EoS_Taylor_AllMu.dat
