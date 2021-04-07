#!/bin/bash


#SBATCH --partition=jn511_1,main             # Partition (job queue)
#SBATCH --job-name=XeIC            # Assign an short name to your job
#SBATCH --array=0-10     #number of subfolders needed filled with 999 events
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=500                   # Real memory (RAM) required (MB)
#SBATCH --time=30:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=logs/Xe_out.out
#SBATCH --error=logs/Xe_out.err
#SBATCH --export=ALL                 # Export you current env to the job env
#SBATCH --constraint=hal # for TRENTO must run on hal, other computers don't work

# set data and working directories
cd /projects/jn511_1/trento/build/src


module load mvapich2/2.1 
module load boost/1.60.0

let i="$SLURM_ARRAY_TASK_ID"

mkdir /projects/jn511_1/trento/build/src/XeXe5440TeVdef_b2/"$i" # needed to make subfolders

srun /projects/jn511_1/trento/build/src/trento -c  Xeb2_25.conf -o XeXe5440TeVdef_b2/"$i"   > runs/Xe_"$i".dat #executes command.  Change .conf and output directory as needed


sleep 3
sacct --format NTasks,MaxRSS,Elapsed,AveRSS,AveCPU -j $SLURM_JOBID
sleep 2
