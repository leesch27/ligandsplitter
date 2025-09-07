"""Provide the primary functions."""
import sys, os
import numpy as np
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from .basefunctions import LigandVariables

vars = LigandVariables()

def get_func_groups(ligs):
    all_func_groups = []
    all_derivs_mols = []
    all_derivs_smiles = []
    for i in ligs:
        mol_groups = []
        mol = Chem.MolFromMol2File(f"data/MOL2_files/{i}_H.mol2", sanitize = False)
        mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol))])))[1]
        mol_renum = Chem.RenumberAtoms(mol, mol_neworder)
        for j in vars.functional_groups:
            k = Chem.MolFromSmarts(j)
            if mol_renum.HasSubstructMatch(k):
                mol_groups.append(1) # 1 corresponds to a functional group being present
                derivs, deriv_smiles = derive(i, j)
                all_derivs_mols.extend(derivs)
                all_derivs_smiles.extend(deriv_smiles)
            else:
                mol_groups.append(0) # 0 corresponds to a functional group being absent
        all_func_groups.append(mol_groups)
    return all_func_groups, all_derivs_mols, all_derivs_smiles

def derive(ligand, group):
    derivative_mols = []
    derivative_smile = []

    mol = Chem.MolFromMol2File(f"data/MOL2_files/{ligand}_H.mol2",sanitize=False)
    old_substr = Chem.MolFromSmarts(group)
    if vars.functional_groups_dict[group] in vars.groups_dict:
        replacement_values = vars.groups_dict[vars.functional_groups_dict[group]]
        for a in replacement_values:
            if a != group:
                new_substr = Chem.MolFromSmiles(a)
                rms = AllChem.ReplaceSubstructs(mol, old_substr, new_substr)
                fragments = Chem.GetMolFrags(rms[0])
                for frag in fragments:
                    ind_frag = Chem.MolFragmentToSmiles(rms[0], frag)
                    if (len(ind_frag) > 25) & (ind_frag != "O[H][H]") & (ind_frag not in derivative_smile):
                        derivative_smile.append(ind_frag)
                        temp_mol = Chem.MolFromSmiles(ind_frag)
                        if temp_mol is not None:
                            derivative_mols.append(temp_mol)
    return derivative_mols, derivative_smile