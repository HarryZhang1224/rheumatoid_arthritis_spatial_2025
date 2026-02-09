#!/bin/bash

# Directory to search for subdirectories
base_dir="/data1/deyk/harry/RA_Xenium/data"

# Find subdirectories containing "RA" and extract the string
subdirs=$(find "$base_dir" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | grep "RA" | sed -n 's/.*\(RA[0-9]*\).*/\1/p')

logdir="/data1/deyk/harry/RA_Xenium/job_logs"
script_dir="/data1/deyk/harry/RA_Xenium/job_scripts/runCLQ"
# Loop through each subdirectory name
for subdir in $subdirs; do
    bsub_command="sbatch --partition=cpu \
        --job-name=CLQ_100_${subdir} \
        --nodes=1 \
        --cpus-per-task=32 \
        --mem-per-cpu=25G \
        --time=05:00:00 \
        --output=$logdir/CLQ_100_${subdir}.log \
        --error=$logdir/CLQ_100_${subdir}.log \
        --wrap=\"bash $script_dir/run.sh -s ${subdir} -r 100\""
    eval $bsub_command
done