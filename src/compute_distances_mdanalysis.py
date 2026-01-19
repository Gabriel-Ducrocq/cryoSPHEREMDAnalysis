import argparse
import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
from MDAnalysis.analysis import distances

parser_arg = argparse.ArgumentParser()
parser_arg.add_argument(
    "--reference_structure",
    type=str,
    required=True,
    help="path to the pdb file containing the structure we want to align to.",
)
parser_arg.add_argument(
    "--structures_path",
    type=str,
    required=True,
    help="path to the pdb file containing the structures from which to extract the distances",
)
parser_arg.add_argument(
    "--output_file", type=str, required=False, help="path to the pdb output file"
)


def compute_distances(reference_path, struct_path, output_file):
	"""
	Computes the distance between two parts of the protein
	:param reference_path: str, path to the structure that will serve to define the topology
	:param struct_path: str, path to the folder containing all the structures
	:param output_file: str, path to the filder that will contain the results.
	"""
	all_distances = []
	univ = mda.Universe(reference_path, struct_path)
	atomA = univ.select_atoms("chainID A and resid 321:503")
	atomB = univ.select_atoms("chainID B and resid 321:503")
	for i, ts in tqdm(enumerate(univ.trajectory)):
		centerA = atomA.center_of_mass()
		centerB = atomB.center_of_mass()

		distance = np.sqrt(np.sum((centerA - centerB)**2))
		all_distances.append(distance)

	np.save(output_file, np.array(all_distances))


if __name__ == "__main__":
    args = parser_arg.parse_args()
    reference_path = args.reference_structure
    struct_path = args.structures_path
    output_file = args.output_file
    compute_distances(reference_path, struct_path, output_file)
