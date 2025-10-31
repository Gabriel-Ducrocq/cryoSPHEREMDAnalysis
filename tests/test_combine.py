import os
import tempfile
from pathlib import Path
from src.combine import (
    combine,
    get_files,
    get_indexes,
    sort_files,
    return_pattern,
    read_file,
    write_file,
)


class TestCombine:
    def test_doc(self):
        assert combine.__doc__

    def test_output(self):
        cwd = os.getcwd()
        path_to_struct = f"{cwd}/../data/structures/"
        path_to_file = f"{cwd}/../data/predicted_combined.pdb"
        combine(path_to_struct, path_to_file)
        with open(path_to_file, "r") as f_pred:
            content_pred = f_pred.readlines()
        with open(f"{cwd}/../data/true_combined.pdb", "r") as f:
            content = f.readlines()

        assert content_pred == content


class TestGetFiles:
    def test_doc(self):
        assert get_files.__doc__

    def test_output(self):
        cwd = os.getcwd()
        path_to_struct = f"{cwd}/../data/structures/"
        assert set(get_files(path_to_struct)) == {
            f"{path_to_struct}structure_z_30058.pdb",
            f"{path_to_struct}structure_z_30119.pdb",
            f"{path_to_struct}structure_z_301765.pdb",
        }

    def test_input(self):
        path_to_struct = "non/existing/path/"
        has_thrown = False
        try:
            get_files(path_to_struct)
        except ValueError:
            has_thrown = True
        assert has_thrown

        cwd = os.getcwd()
        path_to_struct = f"{cwd}/../data/empty_dir/"
        has_thrown = False
        try:
            get_files(path_to_struct)
        except Exception:
            has_thrown = True
        assert has_thrown


class TestSortFiles:
    def test_doc(self):
        assert sort_files.__doc__

    def test_input(self):
        has_thrown = False
        try:
            sort_files(1, [1, 2])
        except TypeError:
            has_thrown = True

        assert has_thrown

    def test_output(self):
        list_files = [
            "structure_z_2.pdb",
            "structure_z_1.pdb",
            "structure_z_4.pdb",
            "structure_z_12945.pdb",
        ]
        list_indexes = [2, 1, 4, 12945]
        sorted_list_indexes = (1, 2, 4, 12945)
        sorted_list_files = (
            "structure_z_1.pdb",
            "structure_z_2.pdb",
            "structure_z_4.pdb",
            "structure_z_12945.pdb",
        )
        sorted_output = sort_files(list_files, list_indexes)
        assert sorted_output == (sorted_list_indexes, sorted_list_files)


class TestReturnPattern:
    def test_doc(self):
        return_pattern.__doc__

    def test_wrong_input(self):
        file = 1
        has_thrown = False
        try:
            return_pattern(file)
        except Exception:
            has_thrown = True

        assert has_thrown

        file = None
        has_thrown = False
        try:
            return_pattern(file)
        except Exception:
            has_thrown = True

        assert has_thrown

    def test_wrong_format(self):
        file = "structure_10.pdb"
        has_thrown = False
        try:
            return_pattern(file)
        except Exception:
            has_thrown = True

        assert has_thrown

    def test_output(self):
        assert return_pattern("structure_z_10394.pdb") == 10394
        assert return_pattern("test/folder/structure_z_10394.pdb") == 10394


class TestGetIndexes:
    def test_doc(self):
        assert get_indexes.__doc__

    def test_output(self):
        list_files = [
            "structure_z_2.pdb",
            "structure_z_1.pdb",
            "structure_z_4.pdb",
            "structure_z_12945.pdb",
        ]
        assert get_indexes(list_files) == [2, 1, 4, 12945]

        list_files = [
            "path/to45/../structure_z_2.pdb",
            "path/to/../structure_z_1.pdb",
            "path/to/../structure_z_4.pdb",
            "path/to/../structure_z_12945.pdb",
        ]
        assert get_indexes(list_files) == [2, 1, 4, 12945]

    def test_input(self):
        list_files = [
            -1,
            "structure_z_1.pdb",
            "structure_z_4.pdb",
            "structure_z_12945.pdb",
        ]
        has_thrown = False
        try:
            assert get_indexes(list_files)
        except Exception:
            has_thrown = True
        assert has_thrown


class TestReadFile:
    def test_doc(self):
        assert read_file.__doc__

    def test_output(self):
        with tempfile.TemporaryDirectory() as tmp_path:
            tmp_path_lib = Path(tmp_path)
            pdb1 = tmp_path_lib / "structure_z_30405.pdb"
            with open(pdb1, "w") as f:
                f.write(
                    "ATOM      1  N   ALA A   1      0.000  0.000  0.000\n ATOM      2  N   ALA A   1      0.000  0.000  0.000\n"
                )

            assert read_file(pdb1) == [
                "ATOM      1  N   ALA A   1      0.000  0.000  0.000\n",
                " ATOM      2  N   ALA A   1      0.000  0.000  0.000\n",
            ]


class TestWriteFile:
    def test_doc(self):
        assert write_file.__doc__

    def test_output(self):
        with tempfile.TemporaryDirectory() as tmp_path:
            tmp_path_lib = Path(tmp_path)
            pdb1 = tmp_path_lib / "structure_z_30405.pdb"
            with open(pdb1, "w") as f:
                f.write("ATOM      1  N   ALA A   1      0.000  0.000  0.000\n")

            content = [
                "ATOM      1  N   ALA A   1      0.000  0.000  0.000\n",
                "ATOM      2  C   ALA A   2      0.000  0.000  1.000\n",
            ]
            file_path = f"{tmp_path}/written.pdb"
            write_file(content, file_path, 1, is_last=False)
            with open(file_path, "r") as f:
                content_read = f.readlines()

            expected_content = ["MODEL      1\n"] + content + ["ENDMDL\n"]
            assert content_read == expected_content

            write_file(content, file_path, 2, is_last=True)
            with open(file_path, "r") as f:
                content_read = f.readlines()

            expected_content = (
                ["MODEL      1\n"]
                + content
                + ["ENDMDL\n"]
                + ["MODEL      2\n"]
                + content
                + ["ENDMDL\n", "END\n"]
            )
            assert content_read == expected_content


"""
def time_it():
    combine(
        "/Users/gabdu45/PycharmProjects/VAECryoEM/data/dataset/ICLR/Guillaume/cryoSPHERE_three/test_struct/",
        "../data/predicted_combined_big.pdb",
    )
"""

# cProfile.run('time_it()')
