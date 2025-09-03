#!/bin/bash

# Directory to search for subdirectories
base_dir="/lila/data/deyk/harry/spatial/RA_Xenium/data/xenium_ranger_baysor"

# Find subdirectories containing "RA" and extract the string
subdirs=$(find "$base_dir" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | grep "RA" | sed -n 's/.*\(RA[0-9]*\).*/\1/p')

logdir="/lila/data/deyk/harry/spatial/RA_Xenium/job_logs"
script_dir="/lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysorQC"
# Loop through each subdirectory name
for subdir in $subdirs; do
    if [ $subdir == "RA443" ]; then
        bsub_command="bsub -q cpuqueue -J BaysorQC_$subdir -n 20 -R \"rusage[mem=10]\" -W 6:00 -o $logdir/BaysorQC_$subdir.log -e $logdir/BaysorQC_$subdir.log 'sh $script_dir/runBaysorQC.sh -s $subdir'"
        eval $bsub_command
    fi
done