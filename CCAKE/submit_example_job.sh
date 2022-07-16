#!/bin/bash

#-------------------------------------------------------------------------------
# set up results directory
RESULTS_DIRECTORY=example_results
mkdir -p $RESULTS_DIRECTORY

sbatch <<EOT
#!/bin/bash
#SBATCH --output="$RESULTS_DIRECTORY/job.out"

#-------------------------------------------------------------------------------
# set up input parameters file save for reference
INPUT_PARAMETERS_FILE=$RESULTS_DIRECTORY/Input_Parameters.inp
cp input/Input_Parameters_example.inp \${INPUT_PARAMETERS_FILE}

#-------------------------------------------------------------------------------
# save state of current git branch used to generate these results
GIT_SNAPSHOT_FILE=$RESULTS_DIRECTORY/git_snapshot.txt
echo "SHA:" `git rev-parse HEAD` > \${GIT_SNAPSHOT_FILE}
echo "BRANCH:" `git rev-parse --abbrev-ref HEAD` >> \${GIT_SNAPSHOT_FILE}
echo "FULL:" >> \${GIT_SNAPSHOT_FILE}
git log -1 >> \${GIT_SNAPSHOT_FILE}

#-------------------------------------------------------------------------------
# run the code
echo "Running command: ./ccake" \${INPUT_PARAMETERS_FILE} $RESULTS_DIRECTORY
./ccake \${INPUT_PARAMETERS_FILE} $RESULTS_DIRECTORY

exit 0
EOT
