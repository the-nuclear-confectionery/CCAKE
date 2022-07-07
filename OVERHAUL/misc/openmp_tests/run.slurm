#!/bin/bash
#SBATCH --job-name=OMP_hello
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=develop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

export OMP_NUM_THREADS=8
./hello_openmp