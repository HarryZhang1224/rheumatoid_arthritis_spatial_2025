import pandas as pd
import numpy as np
import sys
import os

sample_ids = os.listdir("/data/deyk/harry/spatial/RA_Xenium/data/baysor/")
for sample in sample_ids:
    if (sample != "output-XETG00392__0037379__RA443__20240727__202032"):
        continue
    # Read in the transcript file
    molecule_file = f"/lila/data/deyk/harry/spatial/RA_Xenium/data/baysor/{sample}/filtered_transcripts.csv"
    data_frame = pd.read_csv(molecule_file)
    # Get the median transcripts per cell based on 10x nucleus expansion
    median_transcript_per_cell = int(data_frame[~(data_frame['cell_id'] == "UNASSIGNED")].groupby("cell_id").size().median())
    logdir = "/lila/data/deyk/harry/spatial/RA_Xenium/job_logs"
    script_dir = "/lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysor"
    output_dir = f"/lila/data/deyk/harry/spatial/RA_Xenium/data/baysor/{sample}"
    # min molecules per segment set to half of median transcripts per cell
    min_molecules_per_segment=int(median_transcript_per_cell/2)
    
    content = f"""
    [data]
    x = "x_location" # Name of the x column in the input data. Default: "x"
    y = "y_location" # Name of the y column in the input data. Default: "y"
    z = "z_location" # Name of the y column in the input data. Default: "z"
    gene = "feature_name" # Name of gene column in the input data. Default: "gene"
    force_2d = false # Ignores z-column in the data if it is provided
    min_molecules_per_gene = 1 # Minimal number of molecules per gene. Default: 1
    exclude_genes = "" # Comma-separated list of genes or regular expressions to ignore during segmentation. Example: 'Blank*,MALAT1'
    min_molecules_per_cell = {median_transcript_per_cell} # Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Default: 3
    min_molecules_per_segment = {min_molecules_per_segment} # Minimal number of molecules in a segmented region, required for this region to be considered as a possible cell. Default: min-molecules-per-cell / 4
    confidence_nn_id = 0 # Number of nearest neighbors to use for confidence estimation. Default: min-molecules-per-cell / 2 + 1

    [segmentation]
    scale = -1.0 # Negative values mean it must be estimated from `min_molecules_per_cell`
    scale_std = "25%" # Standard deviation of scale across cells. Can be either number, which means absolute value of the std, or string ended with "%" to set it relative to scale. Default: "25%"
    estimate_scale_from_centers = true # Use scale estimate from DAPI if provided. Default: true

    n_clusters = 6 # Number of clusters to use for cell type segmentation. Default: 4
    prior_segmentation_confidence = 0.3 # Confidence of the prior segmentation. Default: 0.2
    iters = 500 # Number of iterations for the cell segmentation algorithm. Default: 500
    n_cells_init = 0 # Initial number of cells

    nuclei_genes = "" # Comma-separated list of nuclei-specific genes. If provided, `cyto-genes` has to be set, as well.
    cyto_genes = "" # Comma-separated list of cytoplasm-specific genes. If provided, `nuclei-genes` has to be set, as well.

    unassigned_prior_label = "UNASSIGNED" # Label for unassigned cells in the prior segmentation. Default: "0"

    [plotting]
    gene_composition_neigborhood = 0 # Number of neighbors (i.e. 'k' in k-NN), which is used for gene composition visualization. Larger numbers leads to more global patterns. Default: estimate from min-molecules-per-cell
    min_pixels_per_cell = 15 # Number of pixels per cell of minimal size, used to estimate size of the final plot. For most protocols values around 7-30 give enough visualization quality. Default: 15
    max_plot_size = 3000 # Maximum size of the molecule plot in pixels. Default: 5000
    ncv_method = "ri" # Method for gene vector estimation for NCV dimensionality reduction. Possible values: `ri` (Random Indexing), `dense` (PCA), `sparse` (Sparse PCA). Default: ri
    """
    with open(f'{output_dir}/config.toml', 'w') as file:
        file.write(content)
    # Submit the job to run baysor
    os.system(f"bsub -q cpuqueue -J baysor_{sample} -n 15 -R \"rusage[mem=58]\" -q cpuqueue -W 72:00 -o {logdir}/baysor_{sample}.log -e {logdir}/baysor_{sample}.log  'sh {script_dir}/runBaysor.sh -f {molecule_file} -c {output_dir}/config.toml -o {output_dir} -s {sample}'")