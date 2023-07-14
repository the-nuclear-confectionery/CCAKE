#!/bin/bash

sbatch <<EOT
#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 03:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10250

# the script to run
SCRIPT=../misc/plotting_scripts/plot_density_distributions.py
DATAFILE=$1
OUTDIRECTORY=\$(dirname \${DATAFILE})/plots/density_distributions

# create results directory
mkdir -p \${OUTDIRECTORY}

# run it
python3 \${SCRIPT} \${DATAFILE} \${OUTDIRECTORY}

exit 0
EOT
