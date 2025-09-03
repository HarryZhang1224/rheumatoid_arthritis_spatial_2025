#!/bin/bash
while getopts s: option
do
  case "${option}" in 
     s) sample=${OPTARG};;
  esac
done

module load miniforge3
source ~/.bashrc
conda activate /data1/deyk/harry/conda_envs/ST_env
Rscript /data1/deyk/harry/RA_Xenium/job_scripts/runNicheDE/all_subtype_M25_combined/run.R $sample
conda deactivate