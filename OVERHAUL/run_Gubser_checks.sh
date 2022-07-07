#!/bin/bash

mkdir -p $1

sbatch <<EOT
#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10250
#SBATCH --output="$1/job.out"

INPUT_PARAMETERS_FILE=Input_Parameters_Gubser_checks.inp

RESULTS_DIRECTORY=$1
mkdir -p \${RESULTS_DIRECTORY}

# save state of current git branch used to generate these results
GIT_SNAPSHOT_FILE=\${RESULTS_DIRECTORY}/git_snapshot.txt
echo "SHA:" `git rev-parse HEAD` > \${GIT_SNAPSHOT_FILE}
echo "BRANCH:" `git rev-parse --abbrev-ref HEAD` >> \${GIT_SNAPSHOT_FILE}
echo "FULL:" >> \${GIT_SNAPSHOT_FILE}
git log -1 >> \${GIT_SNAPSHOT_FILE}

# save settings file for reference
cp \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}

./ccake \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}

exit 0
EOT
