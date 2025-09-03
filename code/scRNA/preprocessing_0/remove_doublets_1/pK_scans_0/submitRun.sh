#!/bin/bash

# Directory to search for subdirectories
base_dir="/data1/deyk/harry/RA_Xenium/data/Ian_old_scRNA_07062025/seurat/"

# Find subdirectories containing "RA" and extract the string
subdirs=$(find "$base_dir" -maxdepth 1 -type f -name '*.qs' -printf '%f\n' | sed 's/\.qs$//')
logdir="/data1/deyk/harry/RA_Xenium/job_logs"
script_dir="/data1/deyk/harry/Mantel_Zhang_et_al_2025_reproducibility/code/scRNA/preprocessing_0/remove_doublets_1/pK_scans_0"
# Loop through each subdirectory name
for subdir in $subdirs; do
    bsub_command="sbatch --partition=cpu,deyk,cpu_highmem \
        --job-name=DoubletFinder_${subdir} \
        --nodes=1 \
        --cpus-per-task=15 \
        --mem-per-cpu=15G \
        --time=1-00:00:00 \
        --output=$logdir/DoubletFinder_${subdir}.log \
        --error=$logdir/DoubletFinder_${subdir}.log \
        --wrap=\"bash $script_dir/run.sh -s ${subdir}\""
    eval $bsub_command
done