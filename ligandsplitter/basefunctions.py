"""Provide the primary functions."""
import sys, os
import random
import requests
import rcsbapi
from rcsbapi.search import AttributeQuery, Attr, TextQuery, ChemSimilarityQuery
import ipywidgets as widgets
from ipywidgets import FileUpload, Dropdown, Text, Layout, Label, Box, HBox, Button, Output
from IPython.display import display

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

    # create test file path/dir, return error if exists
    testPath = os.path.join(dataPath, "test_files")
    try:
        os.mkdir(testPath)
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

def test_pull():
    q1 = AttributeQuery(attribute = "rcsb_entry_info.selected_polymer_entity_types", operator = "exact_match", value = "Protein (only)")
    q2 = AttributeQuery(attribute = "rcsb_polymer_entity.formula_weight", operator = "less_or_equal", value = 300)
    q3 = AttributeQuery(attribute = "pdbx_database_status.pdb_format_compatible", operator = "exact_match", value = "Y")
    query = q1 & q2 & q3
    global result_random
    result_random = list(query())

def on_button_clicked(b):
    global new_name
    with output:
        print("Loading...")
        new_name = random.choice(result_random)
        print(f"PDB ID retieved: {new_name}")

def pull_random():
    lucky = Button(description="I'm feeling lucky",
               disabled=False,
               button_style='',
               tooltip='Click me to generate a random PDB ID',
               icon='check')
    global output
    output = widgets.Output()
    display(lucky, output)
    lucky.on_click(on_button_clicked)

def get_new_name():
        return new_name

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    create_folders()
