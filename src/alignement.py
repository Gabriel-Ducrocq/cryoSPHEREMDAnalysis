import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj


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
    help="path to the pdb file containing the structures to align",
)
parser_arg.add_argument(
    "--output_file", type=str, required=False, help="path to the pdb output file"
)


def align(reference_path, struct_path, output_file):
    """
    Aligns the structures on the C1 atoms.
    :param reference_path: str, path to the pdb file of the structure serving as a reference for the alignement.
    :param struct_path: str, path of the pdb file containing all the structures.
    :param output_file: str, path to the pdb file we want to save the alignement in.
    :return:
    """
    univ = mda.Universe(reference_path, struct_path)
    univ_ref = mda.Universe(reference_path, reference_path)
    alignement = AlignTraj(univ, univ_ref, select="name C1'", filename=output_file)
    alignement.run()


if __name__ == "__main__":
    args = parser_arg.parse_args()
    reference_path = args.reference_structure
    struct_path = args.structures_path
    output_file = args.output_file
    align(reference_path, struct_path, output_file)
