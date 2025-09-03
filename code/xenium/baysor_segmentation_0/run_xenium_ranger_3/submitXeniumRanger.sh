#!/bin/bash

# List directories containing "XETG" and loop through them
logdir="/lila/data/deyk/harry/spatial/RA_Xenium/job_logs"
script_dir="/lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysor"
for dir in $(find /lila/data/deyk/harry/spatial/RA_Xenium/data/baysor -maxdepth 1 -type d -name "*XETG*"); do
  sample=$(basename "$dir")
  sample_id=$(echo "$sample" | grep -o 'RA[0-9]\+')
  dir_to_check="/lila/data/deyk/harry/spatial/RA_Xenium/data/xenium_ranger_baysor/$sample_id/outs"
  file_to_check="/lila/data/deyk/harry/spatial/RA_Xenium/data/baysor/$sample/segmentation_flight_ready.csv"
  # check if xenium range outputs already exist
  if [ ! -d "$dir_to_check" ] && [ -f "$file_to_check" ]; then
    bsub_command="bsub -q cpuqueue -J Xenium_Ranger_$sample -n 20 -R \"rusage[mem=15G]\" -W 10:00 -o $logdir/Xenium_Ranger_$sample.log -e $logdir/Xenium_Ranger_$sample.log 'sh $script_dir/XeniumRanger.sh -s $sample'"
    eval $bsub_command
  fi
done