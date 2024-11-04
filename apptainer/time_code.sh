#!/bin/bash

export OUTPUT_DIR_BASE="../build/output_visc_gubser/serial/sim"
export AVG_TIME_FILE="serial_time.txt"
export TIME_FILE="time.dat"

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

export OUTPUT_DIR_BASE="../build/output_visc_gubser/cpu/sim"
export AVG_TIME_FILE="cpu_time.txt"
export TIME_FILE="time.dat"

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

export OUTPUT_DIR_BASE="../build/output_visc_gubser/gpu/sim"
export AVG_TIME_FILE="gpu_time.txt"
export TIME_FILE="time.dat"

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

