import os
import tempfile
import argparse
import numpy as np
from pathlib import Path
import MDAnalysis as mda


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

    return []

def clustering():
    """
    Function that clusters the structures based on the C_alpha only (excluding the C1)
    :return:
    """

assert clustering.__doc__
assert get_calpha.__doc__
has_thrown = False
try:
    get_calpha("path/to/file.txt")
except:
    has_thrown = True

assert has_thrown


has_thrown = False
try:
    get_calpha(1)
except:
    has_thrown = True

assert has_thrown


has_thrown = False
try:
    get_calpha(None)
except:
    has_thrown = True

assert has_thrown


has_thrown = False
try:
    get_calpha("fake/path/file.pdb")
except:
    has_thrown = True

assert has_thrown


def test_get_c_alpha_output():
    with tempfile.TemporaryDirectory() as tmp_path:
        tmp_path_lib = Path(tmp_path)
        pdb1 = tmp_path_lib / "test.pdb"
        with open(pdb1, "w") as f:
            f.writelines(["MODEL      1\n",
                          "ATOM      1  CA  MET K   1    -105.733  40.745   5.921  1.00  0.00\n",
                          "ATOM      2  CA  GLU K   2    0.000   0.000   0.000   1.00   0.00\n",
                          "ENDMDL\n",
                          "MODEL      2\n",
                          "ATOM      1  CA  MET K   1    2.000   0.000   0.000   1.00   0.00\n",
                          "ATOM      2  CA  GLU K   2    0.000   0.000   0.000   1.00   0.00\n",
                          "ENDMDL\n",
                          "END"])

        #with open("data/predicted_combined.pdb", "r") as f:
        #    c = f.readlines()

        #print(c[0])
        #print(c[1])
        c = get_calpha(str(pdb1))
        assert get_calpha(str(pdb1)) == np.array([[[-105.733,  40.745, 5.921], [0.0, 0.0, 0.0]], [[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]]])
        #ar = get_calpha("data/predicted_combined.pdb")
        #assert ar == np.array(
        #    [[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]]])




