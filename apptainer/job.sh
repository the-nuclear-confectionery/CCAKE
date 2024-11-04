#!/bin/bash
#SBATCH --account=bbkr-delta-cpu
#SBATCH --job-name="visc_gubser"
##SBATCH --output="a.out.%j.%N.out"
#SBATCH --partition=cpu
##SBATCH --mem=208G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  # could be 1 for py-torch
#SBATCH --cpus-per-task=4   # spread out to use 1 core per numa, set to 64 if tasks is 1
##SBATCH --constraint="scratch"
##SBATCH --gpus-per-node=4
##SBATCH --gpu-bind=none   # select a cpu close to gpu on pci bus topology
##SBATCH --account=bbkr-delta-gpu    # <- match to a "Project" returned by the "accounts" command
##SBATCH --exclusive  # dedicated node for this job
#SBATCH --no-requeue
#SBATCH -t 04:00:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

export OMP_NUM_THREADS=4  # if code is not multithreaded, otherwise set to 8 or 16
srun make run-long


