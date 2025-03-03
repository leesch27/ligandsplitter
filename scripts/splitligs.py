import numpy as np
#import pandas as pd
import re
import sys, os
import importlib.util
from Bio.PDB import PDBList, PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
import MDAnalysis as mda 
from openbabel import pybel

#from splitligfunctions import create_folders, retrieve_pdb_file, separate_mol2_ligs
def load_module(name, path):
    dirname = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(dirname, path)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module

splitligfunctions = load_module("splitligfunctions", "splitligfunctions.py")

def main(pdb_id = None, format = None):
    splitligfunctions.create_folders()
    splitligfunctions.retrieve_pdb_file(pdb_id, format = format)
    splitligfunctions.separate_mol2_ligs(filename = f"data/MOL2_files/{pdb_id}_ligand.mol2")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
        sys.exit(1)
    main(pdb_id = sys.argv[1], format = sys.argv[2])