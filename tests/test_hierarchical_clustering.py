import os
import tempfile
import numpy as np
from pathlib import Path
from src.hierarchical_clustering import write_file, clustering, get_calpha


class TestClustering:
    def test_doc(self):
        assert clustering.__doc__

    def test_input(self):
        has_thrown = False
        try:
            clustering("path/to/file.pdb", "path/to/output.npy")
        except ValueError:
            has_thrown = True

        assert has_thrown

    def test_output(self):
        with tempfile.TemporaryDirectory() as tmp_path:
            tmp_path_lib = Path(tmp_path)
            pdb1 = tmp_path_lib / "test.pdb"
            output_file = tmp_path_lib / "output.npy"
            with open(pdb1, "w") as f:
                f.writelines(
                    [
                        "MODEL      1\n",
                        "ATOM      1  C1  MET K   1    -105.733  40.745   5.921  1.00  0.00\n",
                        "ATOM      2  CA  GLU K   2    0.000   0.000   0.000   1.00   0.00\n",
                        "ENDMDL\n",
                        "MODEL      2\n",
                        "ATOM      1  C1  MET K   1    2.000   0.000   0.000   1.00   0.00\n",
                        "ATOM      2  CA  GLU K   2    2.000   0.000   0.000   1.00   0.00\n",
                        "ENDMDL\n",
                        "END",
                    ]
                )

            distance = np.sqrt(
                np.sum(
                    np.sum(
                        (np.array([0.0, 0.0, 0.0]) - np.array([2.000, 0.000, 0.000]))
                        ** 2
                    )
                )
            )
            linkage_mat = np.zeros((1, 4))
            linkage_mat[0, 0] = 0
            linkage_mat[0, 1] = 1
            linkage_mat[0, 2] = distance
            linkage_mat[0, 3] = 2
            clustering(str(pdb1), str(output_file))
            pred_linkage_mat = np.load(str(output_file))
            assert os.path.isfile(str(output_file))
            assert np.all(pred_linkage_mat == linkage_mat)


class TestGetCalpha:
    def test_doc(self):
        assert get_calpha.__doc__

    def test_input(self):
        def scaffold(path):
            has_thrown = False
            try:
                get_calpha(path)
            except (TypeError, ValueError):
                has_thrown = True

            assert has_thrown

        scaffold("path/to/file.txt")
        scaffold("fake/path/file.pdb")
        scaffold(1)
        scaffold(None)

    def test_output(self):
        with tempfile.TemporaryDirectory() as tmp_path:
            tmp_path_lib = Path(tmp_path)
            pdb1 = tmp_path_lib / "test.pdb"
            with open(pdb1, "w") as f:
                f.writelines(
                    [
                        "MODEL      1\n",
                        "ATOM      1  C1  MET K   1    -105.733  40.745   5.921  1.00  0.00\n",
                        "ATOM      2  CA  GLU K   2    0.000   0.000   0.000   1.00   0.00\n",
                        "ENDMDL\n",
                        "MODEL      2\n",
                        "ATOM      1  C1  MET K   1    2.000   0.000   0.000   1.00   0.00\n",
                        "ATOM      2  CA  GLU K   2    2.000   0.000   0.000   1.00   0.00\n",
                        "ENDMDL\n",
                        "END",
                    ]
                )

            positions = get_calpha(str(pdb1))
            assert np.all(positions.shape == (2, 1, 3))
            assert np.all(positions == np.array([[[0.0, 0.0, 0.0]], [[2.0, 0.0, 0.0]]]))


class TestWrite:
    def test_doc(self):
        assert write_file.__doc__

    def test_output(self):
        with tempfile.TemporaryDirectory() as tmp_path:
            tmp_path_lib = Path(tmp_path)
            output_file = tmp_path_lib / "output.npy"
            linkage = np.random.normal(size=(10, 4))
            assert write_file(str(output_file), linkage) is None
            assert os.path.isfile(str(output_file))
            assert np.all(np.load(str(output_file)) == linkage)
