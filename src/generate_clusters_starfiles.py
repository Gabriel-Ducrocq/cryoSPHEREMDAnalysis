import starfile
import argparse
import pandas as pd
import numpy as np



parser_arg = argparse.ArgumentParser()
parser_arg.add_argument(
    "--centroid_distances_file",
    type=str,
    required=True,
    help="path to the file containing the distance from the structures to the centroids",
)

parser_arg.add_argument(
    "--starfile",
    type=str,
    required=True,
    help="path to the file starfile",
)

parser_arg.add_argument(
    "--output_path",
    type=str,
    required=True,
    help="Path to the folder where we save all the starfiles",
)

def save_starfiles(distances_file, output_path, starfile):
    """
    Extract the lines in the starfile corresponding to each cluster
    :param distances_file: str, path to the file containing the distances between the structures and all the centroids
    :param output_path: str, path to folder containing the multiple starfiles we will save.
    :param starfile: str, path to the starfile from which we want to extract the lines
    :return: None
    """
    dist = np.load(distances_file)
    star = starfile.read(starfile)
    clusters = np.argmin(dist, axis=-1)
    all_unique = np.unique(clusters)
    for cluster_number in all_unique:
        df = star["particles"][clusters == cluster_number]
        new_star = {"optics":star["optics"], "particles":df}
        starfile.write(new_star, output_path + f"starfile_cluster_{cluster_number}.star")




if __name__ == "__main__":
    args = parser_arg.parse_args()
    centroid_distances_file = args.centroid_distances_file
    starfile_path = args.starfile
    output_path = args.output_path
    save_starfiles(centroid_distances_file, output_path, starfile_path)