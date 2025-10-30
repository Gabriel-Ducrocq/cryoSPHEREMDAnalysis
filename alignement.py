import glob
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj


parser_arg = argparse.ArgumentParser()
parser_arg.add_argument(
    "--structures_path",
    type=str,
    required=True,
    help="path to the folder containing the structures to cluster",
)
parser_arg.add_argument(
    "--output_file", type=str, required=False, help="path to the pdb output file"
)


def align(struct_path, output_file):
    """
    Aligns the structures on the C1 atoms.
    :param struct_path: str, path of the folder containing all the structures
    :return:
    """
    traj = glob.glob(f"{struct_path}/*pdb")
    univ = mda.Universe(traj[0], traj)
    univ_ref = mda.Universe(traj[0], traj[0])
    alignement = AlignTraj(univ, univ_ref, select="name C1'", filename=output_file)
    alignement.run()


if __name__ == "__main__":
    args = parser_arg.parse_args()
    struct_path = args.structures_path
    output_file = args.output_file
    align(struct_path, output_file)
