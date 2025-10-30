import re


def return_pattern(filepath):
    """
    Returns the index number of the file name.
    :param filepath: str, path of the file
    :return: integer, the index of the file.
    """
    if not isinstance(filepath, str):
        raise TypeError("The file path must of type str")
    if not re.match("structure_z_[0-9]+\\.pdb", filepath):
        raise ValueError("The file names must match the format 'structure_z_[0-9]+\\.pdb'")
    return int(re.search("[0-9]+", filepath).group(0))


def get_indexes(list_file):
    """
    Get the integer part of the filenames, where the filenames are in the format:
    structure_z_i.pdb where i in the integer.
    :param list_file: list of string, path of each file.
    :return: list of integer, containing the integer of each filename in the same order.
    """
    index_list = [return_pattern(f) for f in list_file]
    return index_list


def combine_files(sorted_list):
    """
    Reads each file in the sorted list and combine them in a single pdb
    :param sorted_list: list of string, the path to each file.
    :return: None
    """
    pass


def get_files(path_to_structures):
    """
    Get the list of the files in the structure folder and return it
    :param path_to_structures: str, path to the folder containing the structures.
    :return: list of the files path
    """
    return []

def sort_files(list_file):
    """
    Sort the list of files by the number of the structure. The file names are in the format:
    structure_z_i.pdb where i is an integer.
    :param list_file: list of str, list of path to the pdb.
    :return: list of str, list of path to the pdb but sorted by their structure number.
    """
    return []

def combine(path_to_structures):
    """
    This function combines different pdb files in a single one
    containing different models.
    Arguments:
        - path_to_structures: str, path to the folder containing the structures
    :return:
    """
    list_files = get_files(path_to_structures)
    sorted_list_files = sort_files(list_files)




assert combine.__doc__
assert combine("path/to/structrues/") is None
assert get_files.__doc__
assert sort_files.__doc__
assert combine_files.__doc__
assert get_indexes.__doc__
list_files = ["structure_z_2.pdb", "structure_z_1.pdb", "structure_z_4.pdb", "structure_z_12945.pdb"]
assert get_indexes(list_files) == [2, 1, 4, 12945]
list_files = [-1, "structure_z_1.pdb", "structure_z_4.pdb", "structure_z_12945.pdb"]
has_thrown = False
try:
    assert get_indexes(list_files)
except:
    has_thrown = True
assert has_thrown

assert return_pattern.__doc__



#assert sort_files(list_files) == ["structure_z_1.pdb", "structure_z_2.pdb", "structure_z_3.pdb", "structure_z_4.pdb"]
