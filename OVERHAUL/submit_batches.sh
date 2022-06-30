#!/bin/bash

mkdir -p $1

sbatch <<EOT
#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10250
#SBATCH --array=0-$3
#SBATCH --output="$1/job_%a.out"

echo This is thread = "\${SLURM_ARRAY_TASK_ID}"

cen=$1
mkdir -p \${cen}

INPUT_PARAMETERS_FILE=\${cen}/Input_Parameters_\${SLURM_ARRAY_TASK_ID}.inp

SLURM_SCRIPT_IC_PATH=$2/densities\${SLURM_ARRAY_TASK_ID}.dat

echo "Check: SLURM_SCRIPT_IC_PATH="\${SLURM_SCRIPT_IC_PATH}

cat Input_Parameters_generator.inp \
  | sed "s:SLURM_SCRIPT_IC_PATH:\${SLURM_SCRIPT_IC_PATH}:g" \
  > \${INPUT_PARAMETERS_FILE}

RESULTS_DIRECTORY=\${cen}/results-\${SLURM_ARRAY_TASK_ID}

mkdir -p \${RESULTS_DIRECTORY}

# save state of current git branch used to generate these results
GIT_SNAPSHOT_FILE=\${RESULTS_DIRECTORY}/git_snapshot.txt
echo "SHA:" `git rev-parse HEAD` > \${GIT_SNAPSHOT_FILE}
echo "BRANCH:" `git rev-parse --abbrev-ref HEAD` >> \${GIT_SNAPSHOT_FILE}
echo "FULL:" >> \${GIT_SNAPSHOT_FILE}
git log -1 >> \${GIT_SNAPSHOT_FILE}

# save settings file for reference
cp \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}

./persephone \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}
#echo "Next step: run ./persephone" \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}

exit 0
EOT
