#!/bin/bash                                                                                                                                                             

#1- DIRECTORY
#2- FILENAME
#3- EVENTS
#4- GRID


sed 's,DIRECTORY,'"$1"',; s,FILENAME,'"$1"', ; s,EVENTS,0-'"$3"', ; s,GRID,'"$4"','  ECC.sh > ECC_"$1".sh

chmod +x ECC_"$1".sh

sbatch ECC_"$1".sh
