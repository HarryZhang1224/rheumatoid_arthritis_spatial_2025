#!/bin/bash

# Directory to search for subdirectories
base_dir="/data1/deyk/harry/RA_Xenium/data"

# Find subdirectories containing "RA" and extract the string
subdirs=$(find "$base_dir" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | grep "RA" | sed -n 's/.*\(RA[0-9]*\).*/\1/p')
# subdirs="RA331"
logdir="/data1/deyk/harry/RA_Xenium/job_logs"
script_dir="/data1/deyk/harry/RA_Xenium/job_scripts/runNicheDE/all_subtype_M25_combined"
# Loop through each subdirectory name
for subdir in $subdirs; do
    bsub_command="sbatch --partition=cpu \
        --job-name=nicheDE_all_subtype_M25_combined_${subdir} \
        --nodes=1 \
        --cpus-per-task=10 \
        --mem-per-cpu=50G \
        --time=5-00:00:00 \
        --output=$logdir/nicheDE_all_subtype_M25_combined_${subdir}.log \
        --error=$logdir/nicheDE_all_subtype_M25_combined_${subdir}.log \
        --wrap=\"bash $script_dir/run.sh -s ${subdir}\""
    eval $bsub_command
done