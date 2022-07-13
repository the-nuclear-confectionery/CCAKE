#!/bin/bash

mkdir -p $1

sbatch <<EOT
#!/bin/bash
#SBATCH --output="$1/job.out"

#-------------------------------------------------------------------------------
# set up results directory and input parameters file
RESULTS_DIRECTORY=$1
mkdir -p \${RESULTS_DIRECTORY}

INPUT_PARAMETERS_FILE=\${RESULTS_DIRECTORY}/Input_Parameters.inp

cp input/Input_Parameters_example.inp \${INPUT_PARAMETERS_FILE}

#-------------------------------------------------------------------------------
# save state of current git branch used to generate these results
GIT_SNAPSHOT_FILE=\${RESULTS_DIRECTORY}/git_snapshot.txt
echo "SHA:" `git rev-parse HEAD` > \${GIT_SNAPSHOT_FILE}
echo "BRANCH:" `git rev-parse --abbrev-ref HEAD` >> \${GIT_SNAPSHOT_FILE}
echo "FULL:" >> \${GIT_SNAPSHOT_FILE}
git log -1 >> \${GIT_SNAPSHOT_FILE}

#-------------------------------------------------------------------------------
# save settings file for reference
cp \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}

#-------------------------------------------------------------------------------
# run the code
echo "Running command: ./ccake" \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}
./ccake \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}

exit 0
EOT
