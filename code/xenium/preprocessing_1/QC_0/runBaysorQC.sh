#!/bin/bash
while getopts s: option
do
  case "${option}" in 
     s) sample=${OPTARG};;
  esac
done
rm /lila/data/deyk/harry/spatial/RA_Xenium/job_logs/BaysorQC_$sample.log

source activate /lila/data/deyk/harry/spatial/conda_envs/ST_env
Rscript /lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysorQC/runBaysorQC.R $sample
conda deactivate