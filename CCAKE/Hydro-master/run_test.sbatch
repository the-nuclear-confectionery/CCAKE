#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 7-00:00:00
#SBATCH --mem=128G

mkdir -p outputfiles

/usr/bin/time -v  ./vusphydro settings.inp 0 0
