#!/bin/bash
while getopts s: option
do
  case "${option}" in 
     s) sample=${OPTARG};;
  esac
done

module load miniforge3
source ~/.bashrc
conda activate /data1/deyk/harry/conda_envs/DoubletFinder
Rscript /data1/deyk/harry/Mantel_Zhang_et_al_2025_reproducibility/code/scRNA/preprocessing_0/remove_doublets_1/detect_doublets_1/run.R $sample
conda deactivate