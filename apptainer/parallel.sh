#!/bin/bash

#SBATCH --account=bbkr-delta-gpu
#SBATCH --job-name="vg_para"
##SBATCH --output="a.out.%j.%N.out"
#SBATCH --partition=gpuA100x4
#SBATCH --mem=208G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=closest
#SBATCH --exclusive
#SBATCH -t 30:00:00
#SBATCH --constraint="scratch"
##SBATCH -e slurm-%j.err
##SBATCH -o slurm-%j.out


# Export variables for use in the task script
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=true
export OUTPUT_DIR_BASE="../build/output_visc_gubser/para/sim"
export AVG_TIME_FILE="para_time.log"
export TIME_FILE="time.dat"

# Set variables for simulation
NUM_RUNS=10

# Set variables for simulation
cat << EOF > run_config.sh
OUTPUT_DIR_BASE="$OUTPUT_DIR_BASE"
TIME_FILE="$TIME_FILE"
OMP_NUM_THREADS=$OMP_NUM_THREADS
EOF

# Record start time
start_time=$(date +%s)

# Run the tasks and time it
echo "Starting parallel tasks at $(date)"
{ time seq 1 $NUM_RUNS | xargs -P 64 -I {} ./tasks.sh {} $OUTPUT_DIR_BASE $TIME_FILE ; } 2> time_taken.log

wait

# Record end time
end_time=$(date +%s)
runtime=$((end_time - start_time))
avg_time_per_run=$((runtime / NUM_RUNS))
echo "Total runtime (from bash): $runtime seconds" >> $AVG_TIME_FILE
echo "Average time per run (from bash): $avg_time_per_run seconds" >> $AVG_TIME_FILE

# Calculate averages
echo "Number of cpus per task: $SLURM_CPUS_PER_TASK" >> $AVG_TIME_FILE
echo "Number of gpus per task: $SLURM_GPUS_PER_TASK" >> $AVG_TIME_FILE
echo "Number of tasks per node: $SLURM_NTASKS_PER_NODE" >> $AVG_TIME_FILE
awk '
    { time+=$1; steps+=$2; count++ }
    END {
        avg_time = time/count;
        avg_steps = steps/count;
        avg_time_per_step = avg_time/avg_steps;
        print "Average simulation time: " avg_time " seconds";
        print "Average timesteps per run: " avg_steps;
        print "Average time per timestep: " avg_time_per_step " seconds";
        print "Number of runs: " count;
    }
' $OUTPUT_DIR_BASE*/time.dat >> $AVG_TIME_FILE

echo "All simulations completed. Results stored in $AVG_TIME_FILE."

rm run_config.sh
