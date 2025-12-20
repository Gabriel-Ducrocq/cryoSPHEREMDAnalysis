import utils
import sklearn
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap


parser_arg = argparse.ArgumentParser()
parser_arg.add_argument(
    "--structures_path",
    type=str,
    required=True,
    help="path to the pdb file containing the structure we want to align to",
)
parser_arg.add_argument(
    "--distances_mat",
    type=str,
    required=True,
    help="Path the output npy file containing the linkage matrix",
)


def compute_pairwise_distances(structures_path, output_path):
    """
    Function computing the pairwise distances between structures
    :param structures_path: str, path to the pdb file containing all the structures
    :param output_path: str, path where we want to save the pairwise distances matrix
    :return: None
    """
    univ = mda.Universe(structures_path)
    result = diffusionmap.DistanceMatrix(univ, select='name CA').run()
    utils.write_file(output_path, result.dist_matrix)


if __name__ == "__main__":
    args = parser_arg.parse_args()
    structures_path = args.structures_path
    distances_mat_path = args.distances_mat
    compute_pairwise_distances(structures_path, distances_mat_path)