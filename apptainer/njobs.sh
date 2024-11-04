#!/bin/bash
#SBATCH --account=bbkr-delta-cpu
#SBATCH --job-name="serial"
#SBATCH --partition=cpu
#SBATCH --mem=208G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --constraint="scratch"
##SBATCH --gpus-per-node=4
##SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --no-requeue
export OMP_NUM_THREADS=32

# Set variables for simulation
NUM_RUNS=50                          # Number of runs
TIME_FILE="time.dat"                 # File where simulation time is logged
AVG_TIME_FILE="serial_time.log"     # File to store the average time
TIME_SUM=0                           # Sum of times for all runs
TIMESTEP_SUM=0                       # Sum of timesteps for all runs
OUTPUT_DIR_BASE="../build/output_visc_gubser/serial"             # Base directory for outputs

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
    make --output=out_serial_${i}.out --error=err_serial_${i}.err run-visc-gubser OUTPUT_DIR=$OUTPUT_DIR
    
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

echo "CPUS per task: $SLURM_CPUS_PER_TASK" > AVG_TIME_FILE
echo "Tasks per node: $SLURM_NTASKS_PER_NODE" > AVG_TIME_FILE
# Calculate and store average time per run
AVG_TIME=$(echo "$TIME_SUM / $NUM_RUNS" | bc -l)
echo "Average simulation time: $AVG_TIME seconds" > $AVG_TIME_FILE

# Calculate and store average timesteps per run
AVG_TIMESTEPS=$(echo "$TIMESTEP_SUM / $NUM_RUNS" | bc -l)
echo "Average timesteps per run: $AVG_TIMESTEPS" >> $AVG_TIME_FILE

# Calculate and store average time per timestep
AVG_TIME_PER_TIMESTEP=$(echo "$AVG_TIME / $AVG_TIMESTEPS" | bc -l)
echo "Average time per timestep: $AVG_TIME_PER_TIMESTEP seconds" >> $AVG_TIME_FILE

echo "All simulations completed. Results stored in $AVG_TIME_FILE."

