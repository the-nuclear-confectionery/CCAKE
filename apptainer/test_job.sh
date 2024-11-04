#!/bin/bash
#SBATCH --account=bbkr-delta-cpu
#SBATCH --job-name="test"
#SBATCH --output="a.out.%j.%N.out"
#SBATCH --partition=cpu
#SBATCH --mem=208G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --constraint="scratch"
##SBATCH --gpus-per-node=4
##SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --no-requeue
#SBATCH -t 04:00:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

export OMP_NUM_THREADS=16

# Set variables for simulation
NUM_RUNS=1                          # Number of runs
TIME_FILE="time.dat"                 # File where simulation time is logged
AVG_TIME_FILE="average_time.log"     # File to store the average time
TIME_SUM=0                           # Sum of times for all runs
TIMESTEP_SUM=0                       # Sum of timesteps for all runs
OUTPUT_DIR_BASE="../build/output_visc_gubser"             # Base directory for outputs

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
OUTPUT_DIR="${OUTPUT_DIR_BASE}"
mkdir -p $OUTPUT_DIR

    # Run the makefile target and the executable
echo "Running simulation..."
make run-visc-gubser OUTPUT_DIR=$OUTPUT_DIR
./a.out > "$OUTPUT_DIR/myjob.out"

    # Assume the time and timesteps are logged in a file and move it into the output folder
    #mv $TIME_FILE "$OUTPUT_DIR/"

    # Extract the time and timesteps, add them to running sums
TIME_AND_STEPS=$(extract_time_and_timesteps $OUTPUT_DIR)
TIME=$(echo $TIME_AND_STEPS | awk '{print $1}')
TIMESTEPS=$(echo $TIME_AND_STEPS | awk '{print $2}')

TIME_SUM=$(echo "$TIME_SUM + $TIME" | bc)
TIMESTEP_SUM=$(echo "$TIMESTEP_SUM + $TIMESTEPS" | bc)	

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
