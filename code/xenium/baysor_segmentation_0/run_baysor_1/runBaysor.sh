#!/bin/bash
while getopts f:c:o:s: option
do
  case "${option}" in 
     f) file=${OPTARG};;
     c) config=${OPTARG};;
     o) output=${OPTARG};;
     s) sample=${OPTARG};;
  esac
done

rm /lila/data/deyk/harry/spatial/RA_Xenium/job_logs/baysor_$sample.log
source /lila/home/zhangh10/.bashrc
# Run Baysor after getting the parameter configurations for each sample
JULIA_NUM_THREADS=15 baysor run -c $config --polygon-format GeometryCollection -o $output -p $file :cell_id