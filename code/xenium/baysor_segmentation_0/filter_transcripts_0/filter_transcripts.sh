#!/bin/bash
while getopts s: option
do
  case "${option}" in 
     s) sample=${OPTARG};;
  esac
done

source activate /lila/data/deyk/harry/spatial/conda_envs/envi
python /lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysor/filter_transcripts.py -folder_name ${sample} -min_qv 20