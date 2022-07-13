#!/bin/bash

sbatch <<EOT
#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 03:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10250

# the script to run
SCRIPT=misc/plotting_scripts/animate_profile_from_HDF.py
DATAFILE=$1
QUANTITY=$2
OUTDIRECTORY=\$(dirname \${DATAFILE})/movies

# create results directory
mkdir -p \${OUTDIRECTORY}

# run it
python3 \${SCRIPT} \${DATAFILE} \
                   \${QUANTITY} \
                   \${OUTDIRECTORY}/\${QUANTITY}"_evolution.mp4"

exit 0
EOT
