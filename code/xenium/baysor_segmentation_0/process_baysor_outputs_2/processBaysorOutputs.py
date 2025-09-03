import json
import geojson
import numpy as np
import pandas as pd
import argparse

def main():
    # Parse input arguments.
    args = parse_args()

    # Load GeoJSON object
    with open(f"/data/deyk/harry/spatial/RA_Xenium/data/baysor/{args.folder_name}/segmentation_polygons_2d.json", "r") as geojson:
        data = json.load(geojson)
    # Get all the cell ids in the GeoJSON object
    all_cell_ids = []
    for entry in data["geometries"]:
        all_cell_ids.append(entry["cell"])
    polygon_unique_cell_ids = set(all_cell_ids)
    # Read the segmentation file
    df = pd.read_csv(f"/data/deyk/harry/spatial/RA_Xenium/data/baysor/{args.folder_name}/segmentation.csv")
    seg_cell_ids = set(df[df["cell"].notna()]["cell"])
    # Subset the segmentation file to those drawn in the GeoJSON object, required for Xenium Ranger to run successfully
    exclude_seg_cell = seg_cell_ids.difference(polygon_unique_cell_ids)
    df = df[~df["cell"].isin(set(exclude_seg_cell))]
    # Save the subsetted segmentation file
    df.to_csv(f"/data/deyk/harry/spatial/RA_Xenium/data/baysor/{args.folder_name}/segmentation_flight_ready.csv", index=None)
    # Convert cell ids in the GeoJSON file to integers, required for Xenium Ranger to run successfully
    def convert_cell_id_to_int(data):
        for entry in data:
            cell_id_str = entry['cell'].split('-')[-1]
            cell_id_int = int(cell_id_str)
            entry['cell'] = cell_id_int
    # Convert cell IDs in the data
    convert_cell_id_to_int(data["geometries"])# Save the GeoJSON object with converted cell ids
    with open(f"/data/deyk/harry/spatial/RA_Xenium/data/baysor/{args.folder_name}/segmentation_polygons_2d_flight_ready.json", 'w') as json_file:
        json.dump(data, json_file)


#--------------------------
# Helper functions

def parse_args():
    """Parses command-line options for main()."""
    summary = 'Subset Baysor segmentation results and rename cell ids in GeoJSON object to ensure Xenium Ranger flight success'

    parser = argparse.ArgumentParser(description=summary)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-folder_name',
                               required = True,
                               help="Folder name in the Baysor output directory")

    try:
        opts = parser.parse_args()
    except:
        sys.exit(0)

    return opts



if __name__ == "__main__":
    main()