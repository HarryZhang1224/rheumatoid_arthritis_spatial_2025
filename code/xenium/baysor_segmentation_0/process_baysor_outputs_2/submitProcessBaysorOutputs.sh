#!/bin/bash

# List directories containing "XETG" and loop through them
logdir="/lila/data/deyk/harry/spatial/RA_Xenium/job_logs"
script_dir="/lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysor"
for dir in $(find /lila/data/deyk/harry/spatial/RA_Xenium/data/baysor -maxdepth 1 -type d -name "*XETG*"); do
  sample=$(basename "$dir")
  # check if output files already exist
  file_to_check_1="$dir/segmentation_flight_ready.csv"
  file_to_check_2="$dir/segmentation.csv"
  if [ ! -f "$file_to_check_1" ] && [ -f "$file_to_check_2" ]; then
    bsub_command="bsub -q cpuqueue -J ranger_preflight_process_$sample -n 10 -R \"rusage[mem=15G]\" -W 2:00 -o $logdir/ranger_preflight_process_$sample.log -e $logdir/ranger_preflight_process_$sample.log 'sh $script_dir/processBaysorOutputs.sh -s $sample'"
    eval $bsub_command
  fi
done