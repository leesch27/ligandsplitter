import numpy as np
#import pandas as pd
import re
import sys, os
import importlib.util
from Bio.PDB import PDBList, PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
import MDAnalysis as mda 
from openbabel import pybel
from ligandsplitter.basefunctions import create_folders
from ligandsplitter.ligandsplit import retrieve_pdb_file, separate_mol2_ligs

def main(pdb_id = None, format = None):
    create_folders()
    retrieve_pdb_file(pdb_id, format = format)
    separate_mol2_ligs(filename = f"data/MOL2_files/{pdb_id}_ligand.mol2")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
        sys.exit(1)
    main(pdb_id = sys.argv[1], format = sys.argv[2])