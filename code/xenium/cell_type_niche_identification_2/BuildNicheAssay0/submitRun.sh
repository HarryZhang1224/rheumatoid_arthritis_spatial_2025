#!/bin/bash

# Directory to search for subdirectories
base_dir="/data1/deyk/harry/RA_Xenium/data"

# Find subdirectories containing "RA" and extract the string
subdirs=$(find "$base_dir" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | grep "RA" | sed -n 's/.*\(RA[0-9]*\).*/\1/p')

logdir="/data1/deyk/harry/RA_Xenium/job_logs"
script_dir="/data1/deyk/harry/RA_Xenium/job_scripts/runBuildNicheAssay"
# Loop through each subdirectory name
subdirs="RA480"
for subdir in $subdirs; do
    bsub_command="sbatch --partition=cpu \
        --job-name=runBuildNicheAssay_${subdir} \
        --nodes=1 \
        --cpus-per-task=10 \
        --mem-per-cpu=40G \
        --time=8:00:00 \
        --output=$logdir/runBuildNicheAssay_${subdir}.log \
        --error=$logdir/runBuildNicheAssay_${subdir}.log \
        --wrap=\"bash $script_dir/run.sh -s ${subdir}\""
    eval $bsub_command
    sleep 5
done