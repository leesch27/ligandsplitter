"""Provide the primary functions."""
import sys, os

def create_folders():
    """
    Create data path/dir and sub folders for each file type, return error if exists

    Returns
    -------
    current_dir : str
        Current directory that is in use.
    """

    global current_dir
    current_dir = os.getcwd()
    dataPath = os.path.join(current_dir, "data")
    try:
        os.mkdir(dataPath)
    except OSError as error:
        print(error)

    # create pdb file path/dir, return error if exists
    pdbPath = os.path.join(dataPath, "PDB_files")
    try:
        os.mkdir(pdbPath)
    except OSError as error:
        print(error)

    # create mol2 file path/dir, return error if exists
    mol2Path = os.path.join(dataPath, "MOL2_files")
    try:
        os.mkdir(mol2Path)
    except OSError as error:
        print(error)

    # create pdbqt file path/dir, return error if exists
    pdbqtPath = os.path.join(dataPath, "PDBQT_files")
    try:
        os.mkdir(pdbqtPath)
    except OSError as error:
        print(error)
    return current_dir

def convert_type(start_type):
    """
    Placeholder function to show example docstring (NumPy format).

    Replace this function and doc string for your own project.

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    try:
        isinstance(int(start_type), int)  
        return True
    except ValueError:
        return False

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    create_folders()
