#!/bin/bash
while getopts s: option
do
  case "${option}" in 
     s) sample=${OPTARG};;
  esac
done

rm /lila/data/deyk/harry/spatial/RA_Xenium/job_logs/ranger_preflight_process_$sample.log
source activate /lila/data/deyk/harry/spatial/conda_envs/envi
python /lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysor/processBaysorOutputs.py -folder_name ${sample}