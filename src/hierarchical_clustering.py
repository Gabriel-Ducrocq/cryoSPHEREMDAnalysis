import os
import scipy
import argparse
import numpy as np
import MDAnalysis as mda

parser_arg = argparse.ArgumentParser()
parser_arg.add_argument(
    "--structures_path",
    type=str,
    required=True,
    help="path to the pdb file containing the structure we want to align to",
)
parser_arg.add_argument(
    "--output_path",
    type=str,
    required=True,
    help="Path the output npy file containing the linkage matrix",
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


def clustering(file_path, output_path):
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


if __name__ == "__main__":
    args = parser_arg.parse_args()
    structures_path = args.structures_path
    output_path = args.output_path
    clustering(structures_path, output_path)
