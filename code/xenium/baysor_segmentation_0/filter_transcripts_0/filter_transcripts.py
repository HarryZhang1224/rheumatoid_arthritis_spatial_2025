import pandas as pd
import numpy as np
import argparse
import sys
import os

def main():
    # Parse input arguments.
    args = parse_args()

    data_frame = pd.read_csv(f"/lila/data/deyk/harry/spatial/RA_Xenium/data/{args.folder_name}/transcripts.csv.gz")

    # Filter transcripts. Ignore negative controls
    filtered_frame = data_frame[(data_frame["qv"] >= args.min_qv) &
                                (~data_frame["feature_name"].str.startswith("NegControlProbe_")) &
                                (~data_frame["feature_name"].str.startswith("antisense_")) &
                                (~data_frame["feature_name"].str.startswith("NegControlCodeword_")) &
                                (~data_frame["feature_name"].str.startswith("BLANK_"))]
    # Output filtered transcripts to CSV
    os.system(f"mkdir -p /lila/data/deyk/harry/spatial/RA_Xenium/data/baysor/{args.folder_name}")
    filtered_frame.to_csv(f"/lila/data/deyk/harry/spatial/RA_Xenium/data/baysor/{args.folder_name}/filtered_transcripts.csv",
                          index=False,
                          encoding = 'utf-8')


#--------------------------
# Helper functions

def parse_args():
    """Parses command-line options for main()."""
    summary = 'Filter transcripts from transcripts.csv based on Q-Score threshold \
               Remove negative controls.'

    parser = argparse.ArgumentParser(description=summary)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-folder_name',
                               required = True,
                               help="Folder name in the RA_Xenium data directory with the transcripts.csv.gz file")
    parser.add_argument('-min_qv',
                        default='20.0',
                        type=float,
                        help="The minimum Q-Score to pass filtering. (default: 20.0)")

    try:
        opts = parser.parse_args()
    except:
        sys.exit(0)

    return opts



if __name__ == "__main__":
    main()