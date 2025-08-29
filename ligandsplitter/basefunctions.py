"""Provide the primary functions."""
import sys, os

def create_folders():
    """
    Create data path/dir and sub folders for each file type, return error if exists.

    Parameters
    ----------
    None

    Returns
    -------
    current_dir : String
        Current directory that is in use.
    """

    global current_dir
    current_dir = os.getcwd()
    dataPath = os.path.join(current_dir, "data")
    try:
        os.mkdir(dataPath)
        print("Created directory:", dataPath)
    except OSError as error:
        if error.errno == 17:
            print("Directory exists:", dataPath)
        else:
            print(error)

    # create pdb file path/dir, return error if exists
    pdbPath = os.path.join(dataPath, "PDB_files")
    try:
        os.mkdir(pdbPath)
        print("Created directory:", pdbPath)
    except OSError as error:
        if error.errno == 17:
            print("Directory exists:", pdbPath)
        else:
            print(error)

    # create mol2 file path/dir, return error if exists
    mol2Path = os.path.join(dataPath, "MOL2_files")
    try:
        os.mkdir(mol2Path)
        print("Created directory:", mol2Path)
    except OSError as error:
        if error.errno == 17:
            print("Directory exists:", mol2Path)
        else:
            print(error)

    # create pdbqt file path/dir, return error if exists
    pdbqtPath = os.path.join(dataPath, "PDBQT_files")
    try:
        os.mkdir(pdbqtPath)
        print("Created directory:", pdbqtPath)
    except OSError as error:
        if error.errno == 17:
            print("Directory exists:", pdbqtPath)
        else:
            print(error)

    # create test file path/dir, return error if exists
    testPath = os.path.join(dataPath, "test_files")
    try:
        os.mkdir(testPath)
        print("Created directory:", testPath)
    except OSError as error:
        if error.errno == 17:
            print("Directory exists:", testPath)
        else:
            print(error)
    
    return current_dir

def convert_type(start_type):
    """
    Determine if the start_type is a string or an integer.

    Parameters
    ----------
    start_type : String or int

    Returns
    -------
    bool
        If the start_type is a string, bool is False; if an integer, bool is True.
    """

    try:
        isinstance(int(start_type), int)  
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    create_folders()
