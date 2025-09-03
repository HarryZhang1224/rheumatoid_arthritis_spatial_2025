#!/bin/bash

# Directory to search for subdirectories
base_dir="/data1/deyk/harry/RA_Xenium/data/Ian_old_scRNA_07062025/seurat/"

# Find subdirectories containing "RA" and extract the string
subdirs=$(find "$base_dir" -maxdepth 1 -type f -name '*.qs' -printf '%f\n' | sed 's/\.qs$//')
logdir="/data1/deyk/harry/RA_Xenium/job_logs"
script_dir="/data1/deyk/harry/Mantel_Zhang_et_al_2025_reproducibility/code/scRNA/preprocessing_0/remove_doublets_1/detect_doublets_1"
# Loop through each subdirectory name
for subdir in $subdirs; do
    bsub_command="sbatch --partition=cpushort,deyk,cpu \
        --job-name=DetectDoublets_${subdir} \
        --nodes=1 \
        --cpus-per-task=15 \
        --mem-per-cpu=15G \
        --time=00:30:00 \
        --output=$logdir/DetectDoublets_${subdir}.log \
        --error=$logdir/DetectDoublets_${subdir}.log \
        --wrap=\"bash $script_dir/run.sh -s ${subdir}\""
    eval $bsub_command
done