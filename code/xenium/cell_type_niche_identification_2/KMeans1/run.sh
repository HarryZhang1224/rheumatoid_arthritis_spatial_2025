#!/bin/bash
#SBATCH --job-name="KMeans_niche"
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=40g
#SBATCH --time=12:00:00
#SBATCH --output=/data1/deyk/harry/RA_Xenium/job_logs/KMeans_niche.log
#SBATCH --error=/data1/deyk/harry/RA_Xenium/job_logs/KMeans_niche.log

module load miniforge3
source ~/.bashrc
conda activate /data1/deyk/harry/conda_envs/ST_env
Rscript /data1/deyk/harry/RA_Xenium/job_scripts/runBuildNicheAssay/run.R
conda deactivate