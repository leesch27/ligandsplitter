"""Functions used to retrieve and split ligands from a PDB ID."""
import numpy as np
#import pandas as pd
import re
import sys, os
from Bio.PDB import PDBList, PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
import MDAnalysis as mda 
from openbabel import pybel
from .basefunctions import convert_type
from .ligandgenerate import create_ligands_from_smiles, create_mols_from_smiles
#from ligandsplitter.basefunctions import convert_type
#from ligandsplitter.ligandsmiles import create_ligands_from_smiles, create_mols_from_smiles

class File_Info:
    def __init__(self, tripos_mol, tripos_atom, tripos_bond, lines_mols, lines_atoms, lines_bonds):
        self.tripos_mol = tripos_mol
        self.tripos_atom = tripos_atom
        self.tripos_bond = tripos_bond
        self.lines_mols = lines_mols
        self.lines_atoms = lines_atoms
        self.lines_bonds = lines_bonds

class Ligand:
    def __init__(self, name, lines_atom, lines_bond):
        self.name = name
        self.lines_atom = lines_atom
        self.lines_bond = lines_bond
    @property
    def lines_atom(self):
        return self._lines_atom
    
    @property
    def lines_bond(self):
        return self._lines_bond
        
    @lines_atom.setter
    def lines_atom(self, lines_atom):
        self._lines_atom = lines_atom
        self.num_atoms = lines_atom[-1] - lines_atom[0] + 1
        
    @lines_bond.setter
    def lines_bond(self, lines_bond):
        self._lines_bond = lines_bond
        self.num_bonds = lines_bond[-1] - lines_bond[0] + 1

def retrieve_pdb_file(pdb_id, format = "", type = ""):
    """
    Replace this function and doc string for your own project.

    Parameters
    ----------
    pdb_id : String
        Set to PDB ID of interest
    format: String (pdb or mmcif)

    Returns
    -------
    None
    """

    # isolate protein
    pdb_list = PDBList()
    global pdb_filename
    atom_lines_added = 0
    clean_ligand_exists = True
    # List of residue names for elemental ions
    list_of_ions = ["AG", "AL", "AM", "AU", "AU3", "BA", "BR", "BS3", "CA", "CD", "CE", "CF", "CL", "CO", "3CO", "CR", 
                    "CS", "CU1", "CU", "CU3", "DY", "ER3", "EU3", "EU", "F", "FE", "FE2", "GA", "GD3", "HG", "IN", 'IOD', 
                    "IR3", "IR", "K", "LA", "LI", "LU", "MG", "MN", "MN3", "4MO", "6MO", "NA", "ND", "NI", "3NI", "OS", 
                    "OS4", "PB", "PD", "PR", "PT", "PT4", "4PU", "RB", "RH3", "RHF", "RU", "SB", "SM", "SR", "TB", "TH", 
                    "4TI", "TL", "V", "W", "Y1", "YB", "YB2", "YT3", "ZCM", "ZN", "ZR", "ZTM"]
    if format == "pdb":
        pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir="data/PDB_files", file_format="pdb")
        u = mda.Universe(pdb_filename)
        if type == "protein":
            protein = u.select_atoms("protein")
            protein.write(f"data/PDB_files/{pdb_id}_protein.pdb")
        elif type == "nucleic":
            protein = u.select_atoms("nucleic")
            protein.write(f"data/PDB_files/{pdb_id}_nucleic.pdb")
    
        # isolate ligands and remove water molecules from PDB file
        if type == "protein":
            ligand = u.select_atoms("not protein and not resname HOH")
        elif type == "nucleic":
            ligand = u.select_atoms("not nucleic and not resname HOH")
        try:
            ligand.write(f"data/PDB_files/{pdb_id}_ligand.pdb")
        except IndexError:
            print(f"Protein from PDB ID {pdb_id} has no ligands present. PDB file of protein has been saved to data/PDB_files/{pdb_id}_protein.pdb")
        try:
            with open(f"data/PDB_files/{pdb_id}_clean_ligand.pdb", 'w+') as datafile: 
                with open(f"data/PDB_files/{pdb_id}_ligand.pdb","r") as outfile:
                    data = outfile.readlines()
                for line in data:
                    if 'HETATM' in line:
                        split_line = line.split()
                        line_1_join = split_line[1]
                        line_2_join = split_line[2]
                        line_3_join = split_line[3]
                        if 'HETATM' not in split_line[0]:
                            datafile.write(line)
                        # only write hetatm lines if they are not atomic ions -- if the alphabetical characters in the
                        # res name column and atom name column are the same, it is likely an atomic ion. compare res name
                        # to entries in list_of_ions
                        elif (split_line[0] == 'HETATM') and (line_2_join != line_3_join) and (line_3_join not in list_of_ions):
                            datafile.write(line)
                            atom_lines_added += 1
                        # if res number is 10000 or greater, columns for atom type and res number are counted as one
                        # due to lack of white space, affecting numbering. compare res name to entries in list_of_ions
                        elif (split_line[0] != 'HETATM') and (line_1_join != line_2_join)and (line_2_join not in list_of_ions):
                            datafile.write(line)
                            atom_lines_added += 1
                    else:
                        datafile.write(line)
        except FileNotFoundError:
            clean_ligand_exists = False
            print("")
    elif format == "mmcif":
        pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir="data/PDB_files", file_format="mmCif")
        # parse mmcif file and produce a pdb file with only the protein present
        p = MMCIFParser()
        struc_prot = p.get_structure("", pdb_filename)
        for model in struc_prot:
            for chain in model:
                for residue in list(chain):
                    res_id = residue.id
                    if res_id[0] != ' ':
                        chain.detach_child(res_id)
        io = PDBIO()
        io.set_structure(struc_prot)
        io.save(f"data/PDB_files/{pdb_id}_protein.pdb")
        # get ligand information from mmcif file if available
        p = MMCIFParser()
        struc_lig = p.get_structure("", pdb_filename)
        for model in struc_lig:
            for chain in model:
                for residue in list(chain):
                    res_id = residue.id
                    if res_id[0] == ' ' or res_id[0] == 'W':
                        chain.detach_child(res_id)
                    else:
                        for value in list_of_ions:
                            if res_id[0] == f'H_{value}':
                                chain.detach_child(res_id)
        io = PDBIO()
        io.set_structure(struc_lig)
        io.save(f"data/PDB_files/{pdb_id}_clean_ligand.pdb")
        with open(f"data/PDB_files/{pdb_id}_clean_ligand.pdb","r") as outfile:
            data = outfile.readlines()
            for line in data:
                if 'HETATM' in line:
                    atom_lines_added += 1
        if atom_lines_added == 0:
            clean_ligand_exists = False
    else:
        clean_ligand_exists = False
        print("Invalid format entered. Please enter format as either pdb or mmcif.")
    # convert ligand pdb file to mol2 file, or return a warning if only elemental ions are present
    if atom_lines_added == 0 and clean_ligand_exists:
        print(f"Warning: Ligands cannot be extracted from PDB ID {pdb_id} as only atomic ions are present. PDB file of protein has been saved to data/PDB_files/{pdb_id}_protein.pdb")
    elif clean_ligand_exists:
        pdb_mol2 = [m for m in pybel.readfile(filename = f"data/PDB_files/{pdb_id}_clean_ligand.pdb", format='pdb')][0]
        out_mol2 = pybel.Outputfile(filename = f"data/MOL2_files/{pdb_id}_ligand.mol2", overwrite = True, format='mol2')
        out_mol2.write(pdb_mol2)
        print(f"Comprehensive ligand MOL2 file extracted from PDB ID {pdb_id} has been saved to data/MOL2_files/{pdb_id}_ligand.mol2")
    print("Download completed.")

def get_mol2_info(ligand_file):
    """
    Replace this function and doc string for your own project.

    Parameters
    ----------
    pdb_id : String
        Set whether or not to display who the quote is from.

    Returns
    -------
    None
    """
    tripos_mol = []
    tripos_atom = []
    tripos_bond = []
    lines_mols = []
    lines_atoms = []
    lines_bonds = []
    with open(ligand_file, "r") as outfile:
        data = outfile.readlines()
        for linenum, line in enumerate(data):
            if "@<TRIPOS>MOLECULE" in line:
                tripos_mol.append(linenum)
            if "@<TRIPOS>ATOM" in line:
                tripos_atom.append(linenum)
            if '@<TRIPOS>BOND' in line:
                tripos_bond.append(linenum)
    with open(ligand_file, "r+") as outfile:
        data = outfile.readlines()
        a = 1
        # get lines corresponding to molecule information
        for instance, value in enumerate(tripos_mol):
            for linenum, line in enumerate(data):
                for i in range((linenum >= value) and (linenum < value + 5)):
                    lines_mols.append(line)
        #get lines corresponding to atom information
        for instance, value in enumerate(tripos_atom):
            for linenum, line in enumerate(data):
                for i in range((linenum > value) and (linenum < tripos_bond[instance])):
                    if (convert_type(line.split()[0])) and (len(line.split()) > 7):
                        lines_atoms.append(line)
        # get lines corresponding to bond information
        for instance, value in enumerate(tripos_bond):
            temp_bonds = []
            if len(tripos_atom) > (instance + 1):
                for linenum, line in enumerate(data):
                    for i in range((linenum > value) and (linenum < tripos_atom[instance + 1])):
                        if (convert_type(line.split()[0])) and (len(line.split()) == 4):
                            temp_bonds.append(line)
            else:
                for linenum, line in enumerate(data):
                    for i in range(linenum > value):
                        if (convert_type(line.split()[0])) and (len(line.split()) == 4):
                            temp_bonds.append(line)
            # sort temp_bonds by atom numbers
            temp_bonds.sort(key=lambda x: int(x.split()[1]))
            for num, line in enumerate(temp_bonds):
                temp_bond_split = re.split(r"(\s+)", line)
                temp_bond_split[2] = str(num + 1)
                bond_renumbered = ''.join(str(x) for x in temp_bond_split)
                lines_bonds.append(bond_renumbered)
    return File_Info(tripos_mol, tripos_atom, tripos_bond, lines_mols, lines_atoms, lines_bonds)

def get_ligands(file_info, name_vals = {}):
    """
    Replace this function and doc string for your own project.

    Parameters
    ----------
    pdb_id : String
        Set whether or not to display who the quote is from.
    name_vals : dict
        If being used to split a mol2 file generated from SMILES strings, this dict will
        be used to rename the ligands.

    Returns
    -------
    None
    """
    ligs_temp = []
    lig_loc = []
    atoms = []
    UNL1_count = 0
    # get ligand names and atom locations for each ligand
    for linenum, line in enumerate(file_info.lines_atoms):
        ligand = line
        lig_atom = ligand.split()
        lig1 = str(lig_atom[-2])
        # if a ligand is not in the list of identified ligands and is labeled as 
        # "UNL1", attempt to rename it using the name_vals dictionary
        if 'UNL' in lig1:
            try:
                if (len(name_vals) > 0):
                    keys = list(name_vals.values())
                    lig_atom[-2] = keys[UNL1_count]
                    lig1 = keys[UNL1_count]
                    UNL1_count += 1
            except:
                pass
        # if a ligand is not in the list of identified ligands and is not labeled as 
        # "UNL1", record the line number
        if (lig1 not in ligs_temp) & (lig1 != 'UNL1'):
            ligs_temp.append(lig1)
            lig_loc.append(int(linenum))
        # if the number corresponding to the order of atoms is equal to one, it means
        # these atoms belong to a new ligand, record the line number
        elif (int(lig_atom[0]) == 1):
            ligs_temp.append(lig1)
            lig_loc.append(int(linenum))
        # if a ligand is in the list of identified ligands and is a different ligand 
        # than the one in the line above it, the new ligand is a duplicate of a previously
        # identified ligand, record the line number
        elif ((lig1 in ligs_temp) and (lig1 != ligs_temp[-1])):
            ligs_temp.append(lig1)
            lig_loc.append(int(linenum))
    lig_loc.append(len(file_info.lines_atoms))
    # get list containing the number of atoms present in each ligand
    d = 0
    while d < (len(lig_loc) - 1):
        atoms_1 = lig_loc[d + 1] - lig_loc[d]
        atoms.append(int(atoms_1))
        d += 1
    lig_loc_bond = []
    lig_loc_bond.append(0)
    bonds = []
    # get bond locations for each ligand
    for ligand_number, atom in enumerate(atoms):
        ligand_bonds = []
        for linenum, line in enumerate(file_info.lines_bonds):
            ligand = line
            bond_num = ligand.split()
            bond_atom1 = int(bond_num[1])
            bond_atom2 = int(bond_num[2])
            if (max(bond_atom1, bond_atom2) <= (sum(atoms[:ligand_number]) + atom) and min(bond_atom1, bond_atom2) > sum(atoms[:ligand_number])):
                ligand_bonds.append(linenum)
        if len(ligand_bonds) > 0:
            bonds.append(len(ligand_bonds))
            lig_loc_bond.append(ligand_bonds[-1] + 1)
    # create Ligand instances for each ligand
    ligand_list = []
    for ligand_number, ligand in enumerate(ligs_temp):
        ligand_atom1 = lig_loc[ligand_number]
        ligand_atom2 = lig_loc[ligand_number + 1] - 1
        atoms_for_ligand = [ligand_atom1, ligand_atom2]
        ligand_bond1 = lig_loc_bond[ligand_number]
        ligand_bond2 = lig_loc_bond[ligand_number + 1] - 1
        bonds_for_ligand = [ligand_bond1, ligand_bond2]
        new_lig = Ligand(name = ligand, lines_atom = atoms_for_ligand, lines_bond = bonds_for_ligand)
        ligand_list.append(new_lig)
    return ligand_list

def find_ligands_unique(ligand_list):
    """
    Replace this function and doc string for your own project.

    Parameters
    ----------
    pdb_id : String
        Set whether or not to display who the quote is from.

    Returns
    -------
    None
    """
    # get unique ligands based on ligand names
    ligs_unique = []
    ligs_repeat = []
    for index, templig in enumerate(ligand_list):
        temp_lig_name = templig.name
        if temp_lig_name not in ligs_repeat:
            ligs_unique.append(ligand_list[index])
            ligs_repeat.append(temp_lig_name)
    return ligs_unique

def write_mol2(ligs_unique, file_info):
    """
    Replace this function and doc string for your own project.

    Parameters
    ----------
    pdb_id : String
        Set whether or not to display who the quote is from.

    Returns
    -------
    None
    """
    global ligs # name of ligands
    global filenames # resulting file names for each ligand
    ligs = []
    filenames = []
    previous_atoms = 0
    for unique_ind, unique_lig in enumerate(ligs_unique):
        ligs.append(unique_lig.name)
        filename = "data/MOL2_files/" + str(unique_lig.name) + ".mol2"
        filenames.append(filename)
        infile = open(filename, "w") 
        tripos_mols = []
        # write molecule info for ligand
        for line in file_info.lines_mols:
            tripos_mols.append(line)
        temp_mol = re.split(r"(\s+)", tripos_mols[2])
        temp_mol[2] = unique_lig.num_atoms
        temp_mol[4] = unique_lig.num_bonds
        new_mol = ''.join(str(x) for x in temp_mol)
        tripos_mols[2] = new_mol
        tripos_mols.append("\n")
        # write atoms for ligand
        tripos_atoms = ["@<TRIPOS>ATOM\n"]
        counter_atom = 1
        atom1 = unique_lig.lines_atom[0]
        atom2 = unique_lig.lines_atom[-1]
        while atom1 <= atom2:
            temp_atom = re.split(r"(\s+)", file_info.lines_atoms[atom1])
            temp_atom = temp_atom[:-1]
            original_values = file_info.lines_atoms[atom1].split()
            temp_atom[2] = str(counter_atom)
            if len(temp_atom[2]) > 1 and (len(temp_atom[2]) != len(original_values[0])):
                len_space = len(temp_atom[1])
                temp_atom[1] = temp_atom[1][:(len_space + 1 - len(temp_atom[2]))]
            new_atom = ''.join(str(x) for x in temp_atom)
            tripos_atoms.append(new_atom)
            counter_atom += 1
            atom1 += 1
        # write bonds for ligand
        tripos_bonds = ["@<TRIPOS>BOND\n"]
        counter_bond = 1
        bond1 = unique_lig.lines_bond[0]
        bond2 = unique_lig.lines_bond[1]
        while bond1 <= bond2:
            temp_bond = re.split(r"(\s+)", file_info.lines_bonds[bond1])
            temp_bond = temp_bond[:-1]
            original_values = file_info.lines_bonds[bond1].split()
            temp_bond[2] = str(counter_bond)
            if (len(temp_bond[2]) != len(original_values[0])):
                len_diff = len(original_values[0]) - len(temp_bond[2])
                len_space = len(temp_bond[1])
                if len_diff > 0:
                    temp_bond[1] = temp_bond[1] + (" " * (len_diff))
                else:
                    temp_bond[1] = temp_bond[1][:(len_space + 1 - len_diff)]
            temp_bond[4] = str(int(temp_bond[4]) - previous_atoms)
            if (len(temp_bond[4]) != len(original_values[1])):
                len_diff = len(original_values[1]) - len(temp_bond[4])
                len_space = len(temp_bond[3])
                if len_diff > 0:
                    temp_bond[3] = temp_bond[3] + (" " * (len_diff))
                else:
                    temp_bond[3] = temp_bond[3][:(len_space + 1 - len_diff)]
            temp_bond[6] = str(int(temp_bond[6]) - previous_atoms)
            if (len(temp_bond[6]) != len(original_values[2])):
                len_diff = len(original_values[2]) - len(temp_bond[6])
                len_space = len(temp_bond[5])
                if len_diff > 0:
                    temp_bond[5] = temp_bond[5] + (" " * (len_diff))
                else:
                    temp_bond[5] = temp_bond[5][:(len_space + 1 - len_diff)]
            new_bond = ''.join(str(x) for x in temp_bond)
            tripos_bonds.append(new_bond)
            counter_bond += 1
            bond1 += 1
        previous_atoms = previous_atoms + unique_lig.num_atoms
        # write file
        infile.writelines(tripos_mols)
        infile.writelines(tripos_atoms)
        infile.writelines(tripos_bonds)
        infile.close()
    return ligs, filenames

def separate_mol2_ligs(filename = '', name_vals = {}):
    """
    Replace this function and doc string for your own project.

    Parameters
    ----------
    pdb_id : String
        Set whether or not to display who the quote is from.

    Returns
    -------
    None
    """
    current_dir = os.getcwd()
    ligand_file = os.path.join(current_dir, filename)
    file_info = get_mol2_info(ligand_file)
    ligand_list = get_ligands(file_info, name_vals)
    ligs_unique = find_ligands_unique(ligand_list)
    ligs, filenames = write_mol2(ligs_unique, file_info)
    return ligs, filenames
