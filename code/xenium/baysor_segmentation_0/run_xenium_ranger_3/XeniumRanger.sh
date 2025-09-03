#!/bin/bash
while getopts s: option
do
  case "${option}" in 
     s) sample=${OPTARG};;
  esac
done

rm /lila/data/deyk/harry/spatial/RA_Xenium/job_logs/Xenium_Ranger_$sample.log
sample_id=$(echo "$sample" | grep -o 'RA[0-9]\+')
source ~/.bashrc
cd /lila/data/deyk/harry/spatial/RA_Xenium/data/xenium_ranger_baysor

# Run xeniumranger to import baysor results
xeniumranger import-segmentation --xenium-bundle /lila/data/deyk/harry/spatial/RA_Xenium/data/$sample \
--id $sample_id \
--transcript-assignment /lila/data/deyk/harry/spatial/RA_Xenium/data/baysor/$sample/segmentation_flight_ready.csv \
--viz-polygons /lila/data/deyk/harry/spatial/RA_Xenium/data/baysor/$sample/segmentation_polygons_2d_flight_ready.json \
--units microns \
--jobmode local \
--localcores 20 \
--localmem 300 \



