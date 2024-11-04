TIME_FILE="time.dat"
OUTPUT_DIR_BASE="../build/output_visc_gubser"
extract_time_and_timesteps() {
    # Assumes time and timesteps are on the same line, separated by space
    # Example of line in myjob_time.log: "1.23 1000"
    local time_and_steps=$(head -n 1 "$1/$TIME_FILE")
    local time=$(echo $time_and_steps | awk '{print $1}')
    local timesteps=$(echo $time_and_steps | awk '{print $2}')
    echo "$time $timesteps"
}
OUTPUT_DIR="${OUTPUT_DIR_BASE}_1"
TIME_AND_STEPS=$(extract_time_and_timesteps $OUTPUT_DIR)
echo $TIME_AND_STEPS
