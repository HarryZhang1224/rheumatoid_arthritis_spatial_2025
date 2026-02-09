#!/bin/bash
while getopts s:r: option
do
  case "${option}" in 
     s) sample=${OPTARG};;
     r) radius=${OPTARG};;
  esac
done

module load miniforge3
source ~/.bashrc
conda activate /data1/deyk/harry/conda_envs/ST_env
Rscript /data1/deyk/harry/RA_Xenium/job_scripts/runCLQ/run.R $sample $radius
conda deactivate