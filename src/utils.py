import os
import numpy as np
import MDAnalysis as mda



def write_file(path_to_file, pairwise_path):
    """
    Write the linkage matrix obtained as a result of the Hierarchical clustering
    :param path_to_file: str, path to the file we want to save our results in.
    :param pairwise_path: np.array, matrix to save.
    :return: None
    """
    if not isinstance(path_to_file, str) or not path_to_file.endswith(".npy"):
        raise TypeError(
            f"""The path to the output file must be a string with .npy ending", currently {path_to_file}"""
        )
    if len(pairwise_path.shape) != 2:
        raise ValueError(
            f"The linkage matrix must have 2 dimensions, currently it has {len(pairwise_path.shape)}"
        )

    np.save(path_to_file, pairwise_path)



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