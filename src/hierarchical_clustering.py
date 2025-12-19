import os
import scipy
import pickle
import sklearn
import argparse
import numpy as np
from time import time
import MDAnalysis as mda

parser_arg = argparse.ArgumentParser()
parser_arg.add_argument(
    "--structures_path",
    type=str,
    required=True,
    help="path to the pdb file containing the structure we want to align to",
)
parser_arg.add_argument(
    "--linkage_mat_path",
    type=str,
    required=True,
    help="Path the output npy file containing the linkage matrix",
)

parser_arg.add_argument(
    "--model_path",
    type=str,
    required=True,
    help="Path the output pickle file containing the kmeans if needed",
)

parser_arg.add_argument(
    "--algo",
    type=str,
    required=True,
    help="Clutering algorithm to use. Currently supporting only kmeans and hierarchical",
)

parser_arg.add_argument(
    "--n_clusters",
    type=int,
    required=False,
    help="Number of clusters to use.",
)

parser_arg.add_argument(
    "--n_clusters_min",
    type=int,
    required=False,
    help="Minimum number of clusters",
)

parser_arg.add_argument(
    "--n_clusters_max",
    type=int,
    required=False,
    help="Maximum number of clusters",
)

def write_file(path_to_file, linkage):
    """
    Write the linkage matrix obtained as a result of the Hierarchical clustering
    :param path_to_file: str, path to the file we want to save our results in.
    :param linkage: np.array, matrix to save.
    :return: None
    """
    if not isinstance(path_to_file, str) or not path_to_file.endswith(".npy"):
        raise TypeError(
            f"""The path to the output file must be a string with .npy ending", currently {path_to_file}"""
        )
    if len(linkage.shape) != 2:
        raise ValueError(
            f"The linkage matrix must have 2 dimensions, currently it has {len(linkage.shape)}"
        )

    np.save(path_to_file, linkage)


def get_calpha(file_path):
    """
    Function that reads all the structures from a single pdb file and returns their atom positions as a numpy array.
    :param file_path: str, path to the pdb file containing the structures.
    :return: np.array(N_structures, N_c_alpha, 3)
    """
    if not isinstance(file_path, str):
        raise TypeError("The argument file_path must be a string.")

    if not file_path.endswith(".pdb"):
        raise TypeError("The file must be a .pdb file.")

    if not os.path.isfile(file_path):
        raise ValueError(f"The file {file_path} does not exist")

    univ = mda.Universe(file_path, file_path)
    all_coordinates = []
    for _ in univ.trajectory:
        all_coordinates.append(univ.select_atoms("name CA").positions)

    return np.array(all_coordinates)


def clustering(file_path, linkage_mat_path, model_path, n_clusters, algorithm, n_min_clusters=None, n_max_clusters=None):
    """
    Function that clusters the structures based on the C_alpha only (excluding the C1)
    :param file_path: str, path to the pdb file containing the structures.
    :param linkage_mat_path: str, path to the npy file containing the linkage matrix.
    :param model_path: str, path to the pkl file containing the model.
    :param algorithm: str, name of the algorithm we use.
    :param n_clusters: integer, number of clusters.
    :return: None
    """
    assert algorithm in ["kmeans", "hierarchical"]
    if algorithm == "kmeans":
        clustering_kmeans(file_path, linkage_mat_path, model_path, n_clusters, n_min_clusters, n_max_clusters)
    elif algorithm == "hierarchical":
        clustering_hierarchical(file_path, linkage_mat_path)

def clustering_hierarchical(file_path, output_path):
    """
    Function that clusters the structures based on the C_alpha only (excluding the C1)
    :param file_path: str, path to the pdb file containing the structures.
    :param output_path: str, path to the npy file containing the linkage matrix.
    :return: None
    """
    calphas_positions = get_calpha(file_path)
    n_structures = calphas_positions.shape[0]
    calphas_positions_flattened = calphas_positions.reshape(n_structures, -1)
    linkage_matrix = scipy.cluster.hierarchy.linkage(
        calphas_positions_flattened, method="average"
    )
    write_file(output_path, linkage_matrix)


def clustering_kmeans(file_path, linkage_mat_path, model_path, n_cluster, n_min_clusters, n_max_clusters):
    """
    Function that clusters the structures based on the C_alpha only (excluding the C1) using kmeans.
    :param file_path: str, path to the pdb file containing the structures.
    :param linkage_mat_path: str, path to the npy file containing the linkage matrix.
    :param model_path: str, path to the pickled kmeans object.
    :param n_cluster: integer, number of clusters to use with kmeans.
    :return: None
    """
    calphas_positions = get_calpha(file_path)
    n_structures = calphas_positions.shape[0]
    calphas_positions_flattened = calphas_positions.reshape(n_structures, -1)
    if n_min_clusters is None or n_max_clusters is None:
        kmeans_algo = sklearn.cluster.KMeans(n_cluster, verbose=True)
        all_distances = kmeans_algo.fit_transform(calphas_positions_flattened)
        write_file(linkage_mat_path, all_distances)
        with open(model_path, "wb") as f:
            pickle.dump(kmeans_algo, f)
    else:
        for n_cluster in range(n_min_clusters, n_max_clusters+1):
            start = time()
            kmeans_algo = sklearn.cluster.KMeans(n_cluster, verbose=True)
            all_distances = kmeans_algo.fit_transform(calphas_positions_flattened)
            end = time()
            print(f"Time for iteration {n_cluster}", end-start)
            write_file(linkage_mat_path + f"linkage_mat_{n_cluster}.npy", all_distances)
            with open(model_path + f"kmeans_{n_cluster}.pkl", "wb") as f:
                pickle.dump(kmeans_algo, f)






if __name__ == "__main__":
    args = parser_arg.parse_args()
    structures_path = args.structures_path
    linkage_mat_path = args.linkage_mat_path
    model_path = args.model_path
    algorithm = args.algo
    n_clusters = args.n_clusters
    n_max_clusters = args.n_clusters_max
    n_min_clusters = args.n_clusters_min
    clustering(structures_path, linkage_mat_path, model_path, n_clusters, algorithm, n_min_clusters, n_max_clusters)
