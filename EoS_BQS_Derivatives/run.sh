#!/bin/bash
#SBATCH --job-name=EOSbsq.job
#SBATCH --output=EOSbsq.out
#SBATCH --error=EOSbsq.err
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH -A qgp
#SBATCH --partition=qgp

./EoS_BQS Coefficients_Parameters.dat

