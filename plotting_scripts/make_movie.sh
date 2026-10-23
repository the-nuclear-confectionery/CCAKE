#!/bin/bash

# the script to run
echo "Running from directory:" `pwd`
SCRIPT=misc/plotting_scripts/animate_profile_from_HDF.py
DATAFILE=$1
QUANTITY=$2
OUTDIRECTORY=$(dirname ${DATAFILE})/movies

# create results directory
mkdir -p ${OUTDIRECTORY}

# run it
python3 ${SCRIPT} ${DATAFILE} \
                  ${QUANTITY} \
                  ${OUTDIRECTORY}/${QUANTITY}"_evolution.mp4"

exit 0
