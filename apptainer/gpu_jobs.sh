#!/bin/bash
#SBATCH --account=bbkr-delta-gpu
#SBATCH --job-name="vg_gpu"
##SBATCH --output="a.out.%j.%N.out"
#SBATCH --partition=gpuA100x4
#SBATCH --mem=208G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=closest
#SBATCH --exclusive
#SBATCH -t 30:00:00
#SBATCH --constraint="scratch"
##SBATCH -e slurm-%j.err
##SBATCH -o slurm-%j.out

export OMP_NUM_THREADS=64
export OMP_PROC_BIND=true

# Set variables for simulation
NUM_RUNS=1                          # Number of runs
TIME_FILE="time.dat"                 # File where simulation time is logged
AVG_TIME_FILE="gpu_time.log"     # File to store the average time
TIME_SUM=0                           # Sum of times for all runs
TIMESTEP_SUM=0                       # Sum of timesteps for all runs
OUTPUT_DIR_BASE="../build/output_visc_gubser/gpu/sim"             # Base directory for outputs

# Function to extract time and timesteps from the time log file
extract_time_and_timesteps() {
    # Assumes time and timesteps are on the same line, separated by space
    # Example of line in myjob_time.log: "1.23 1000"
    local time_and_steps=$(head -n 1 "$1/$TIME_FILE")
    local time=$(echo $time_and_steps | awk '{print $1}')
    local timesteps=$(echo $time_and_steps | awk '{print $2}')
    echo "$time $timesteps"
}

# Loop over the number of runs
for i in $(seq 1 $NUM_RUNS); do
    OUTPUT_DIR="${OUTPUT_DIR_BASE}_$i"
    mkdir -p $OUTPUT_DIR

    # Run the makefile target and the executable
    echo "Running simulation $i..."
    srun --output=out_gpu_${i}.out --error=err_gpu_${i}.err make run-visc-gubser OUTPUT_DIR="${OUTPUT_DIR_BASE}_$i"

    # Assume the time and timesteps are logged in a file and move it into the output folder
    #mv $TIME_FILE "$OUTPUT_DIR/"

    # Extract the time and timesteps, add them to running sums
    TIME_AND_STEPS=$(extract_time_and_timesteps $OUTPUT_DIR)
    TIME=$(echo $TIME_AND_STEPS | awk '{print $1}')
    TIMESTEPS=$(echo $TIME_AND_STEPS | awk '{print $2}')

    TIME_SUM=$(echo "$TIME_SUM + $TIME" | bc)
    TIMESTEP_SUM=$(echo "$TIMESTEP_SUM + $TIMESTEPS" | bc)
done

wait
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
