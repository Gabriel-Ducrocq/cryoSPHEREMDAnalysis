import starfile
import argparse
import numpy as np
import pandas as pd
import pickle as pkl



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

parser_arg.add_argument(
    "--ctf_path",
    type=str,
    required=True,
    help="Path to the file containing the ctf",
)

parser_arg.add_argument(
    "--poses_path",
    type=str,
    required=True,
    help="Path to the file containing the poses",
)

def save_ctf_poses_pkl(ctf_array, poses_array, output_path, cluster_number, clusters):
    """
    Extract and saves the poses and ctf files corresponding to this particular cluster
    :param ctf_array: np.array containing all the ctf values
    :param poses_array: np.array containing all the poses values
    :param output_path: str, path to the folder to save the results
    :param cluster_number: int, current cluster number
    :param clusters: np.array(int), cluster number for each structure
    :return:
    """
    ctf_cluster = ctf_array[clusters == cluster_number]
    poses_cluster = poses_array[clusters == cluster_number]
    with open(output_path + f"ctf_{cluster_number}.pkl", "wb") as f:
        pkl.dump(ctf_cluster, f)

    with open(output_path + f"poses_{cluster_number}.pkl", "wb") as f:
        pkl.dump(poses_cluster, f)




def save_starfile(cluster_number, clusters, output_path, star):
    """
    Extract the lines in the starfile corresponding to each cluster
    :param cluster_number: np.array(int), containing the current cluster number.
    :param clusters: np.array(int), containing the cluster label for each image
    :param output_path: str, path to folder containing the multiple starfiles we will save.
    :param star: dict of pandas DataFrames corresponding to the starfile
    :return: None
    """
    df = star["particles"][clusters == cluster_number]
    new_star = {"optics":star["optics"], "particles":df}
    starfile.write(new_star, output_path + f"starfile_cluster_{cluster_number}.star")


def treat_clusters(distances_file, output_path, starfile_path, ctf_path, poses_path):
    """
    For each cluster, extract the corresponding line in the starfile and save them in a separate starfile, extracts
    the correspondin ctf and poses values and saves them in pkl files.
    :param distances_file: str, path to the file containing the distances between the structures and all the centroids
    :param output_path: str, path to folder containing the multiple starfiles we will save.
    :param starfile_path: str, path to the starfile from which we want to extract the lines
    :param ctf_path: str, path to the pkl file containing the ctf values
    :param poses_path: str, path to the pkl file containing the poses values
    :return: None
    """
    dist = np.load(distances_file)
    star = starfile.read(starfile_path)
    clusters = np.argmin(dist, axis=-1)
    all_unique = np.unique(clusters)
    with open(ctf_path, "rb") as f:
        ctf_array = pkl.load(f)

    with open(poses_path, "rb") as f:
        poses_array = pkl.load(f)

    for cluster_number in all_unique:
        save_starfile(cluster_number, clusters, output_path, star)
        save_ctf_poses_pkl(ctf_array, poses_array, output_path, cluster_number, clusters)



if __name__ == "__main__":
    args = parser_arg.parse_args()
    centroid_distances_file = args.centroid_distances_file
    starfile_path = args.starfile
    output_path = args.output_path
    ctf_path = args.ctf_path
    poses_path = args.poses_path
    treat_clusters(centroid_distances_file, output_path, starfile_path)