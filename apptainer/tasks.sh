#!/bin/bash

# This script runs a single task
# Usage: ./run_task.sh <run_number>

# Source the configuration file
source run_config.sh

# Get the run number from the command line argument
RUN_NUMBER=$1

# Use the run number to create a unique output directory for this simulation
OUTPUT_DIR="${OUTPUT_DIR_BASE}_${RUN_NUMBER}"
mkdir -p $OUTPUT_DIR

echo "Running simulation $RUN_NUMBER..."
srun --output=out_para_${RUN_NUMBER}.out --error=err_para_${RUN_NUMBER}.err make run-visc-gubser OUTPUT_DIR="$OUTPUT_DIR"

# Check if the time.dat file was created in the output directory
if [ -f "$OUTPUT_DIR/$TIME_FILE" ]; then
    echo "Simulation $RUN_NUMBER completed. Time data collected in $OUTPUT_DIR/$TIME_FILE"
else
    echo "Warning: $TIME_FILE not found in $OUTPUT_DIR after simulation $RUN_NUMBER"
fi
