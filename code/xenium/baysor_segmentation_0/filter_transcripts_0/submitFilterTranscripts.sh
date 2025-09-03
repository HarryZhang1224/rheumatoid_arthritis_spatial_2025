#!/bin/bash

# List directories containing "XETG" and loop through them
logdir="/lila/data/deyk/harry/spatial/RA_Xenium/job_logs"
script_dir="/lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysor"
# Submit a job to filter transcripts for each sample
for dir in $(find /lila/data/deyk/harry/spatial/RA_Xenium/data/ -maxdepth 1 -type d -name "*XETG*"); do
  sample=$(basename "$dir")
  bsub_command="bsub -q cpuqueue -J filter_baysor_$sample -n 10 -R \"rusage[mem=10G]\" -W 3:00 -o $logdir/filter_baysor_$sample.log -e $logdir/filter_baysor_$sample.log 'sh $script_dir/filter_transcripts.sh -s $sample'"
  eval $bsub_command
done