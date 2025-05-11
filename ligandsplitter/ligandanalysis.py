"""Provide the primary functions."""
import sys, os
from rdkit import Chem

def get_vars():
    """
    Initialize essential variables for ligand analysis

    Returns
    -------
    None
    """
    global all_func_groups
    all_func_groups = [] 

    # Make dictionary of residues and residue type, where 1 = hydrophobic alipathic, 2 = hydrophobic aromatic, 
    # 3 = polar uncharged, 4 = polar acidic, and 5 = polar basic
    global type_dict
    type_dict = {
        'GLY': 1,
        'ALA': 1,
        'VAL': 1,
        'LEU': 1,
        'ILE': 1,
        'PRO': 1,
        'MET': 1,
        'PHE': 2,
        'TYR': 2,
        'TRP': 2,
        'SER': 3,
        'THR': 3,
        'CYS': 3,
        'ASN': 3,
        'GLN': 3,
        'ASP': 4,
        'GLU': 4,
        'ARG': 5,
        'HIS': 5,
        'LYS': 5
    }

    #common functional groups
    ester = '[CX3H1,CX3](=O)[OX2H0][C;!$([C]=[O])]'
    ether = '[OD2]([C;!$([C]=[O])])[C;!$([C]=[O])]'
    hydroxy = '[OX2H;$([O]([C])[H]);!$([O](C=O)[H])][H]'
    carbox_acid = '[CX3;!$([CX3][H])](=O)[OX2H][H]'
    aldehyde = '[CX3;$(C([H])(=O)[C])](=[O;!$([O][O])])[H]'
    anhydr = 'C-[CX3](=O)[O][CX3](=O)-C'

    #nitrogen containing
    amine = '[C][NX3;H2;!$(NC=O)]([H])[H]'
    amine_2 = '[C][NX3;H;!$(NC=O)]([C])[H]'
    amine_3 = '[C][NX3;H0;!$(NC=O);!$(N=O)]([CX4])[CX4]'
    amide = '[CX3](=O)[NX3;!$(N=O)]([H])[H]' 
    amide_2 = '[CX3](=O)[NX3;!$(N=O)]([C])[H]'
    amide_3 = '[CX3](=O)[NX3;!$(N=O)]([C])[C]'
    nitro = 'C-N(=O)-O'
    imine = '[C;$(C([C,H])([C,H])=[N])]=[N][C,H]'

    #halogens
    f_hal = 'F'
    cl_hal = 'Cl'
    br_hal = 'Br'
    i_hal = 'I'

    #multiply bonded and rings
    alkene = 'C=C'
    alkyne = 'C#C'
    alkyne_term = 'C#C-[H]'
    phenyl = '[CX4;$([C][c;$(c1cc[c]cc1)]);!$([C]([H])([H])([C])[c;$(c1cc[c]cc1)])][c]1:[c]:[c]:[c]:[c]:[c]1'
    benzyl = '[CX4;$([C](c1ccccc1)([H])([H])[C])][c]1:[c]:[c]:[c]:[c]:[c]1'
    pyrrole = '[$([NH]1C=CC=C1)]'
    imidiz = '[$([NH]1C=NC=C1)]'
    pyridine = '[$([nR1]1[cR1][cR1][cR1][cR1][cR1]1)]'
    pyrimidine = '[$([nR1]1[cR1][nR1][cR1][cR1][cR1]1)]'

    global functional_groups
    functional_groups = [ester, ether, hydroxy, carbox_acid, aldehyde, anhydr, amine, amine_2, amine_3, 
                        amide, amide_2, amide_3, nitro, imine, f_hal, cl_hal, br_hal, i_hal, alkene, 
                        alkyne, alkyne_term, phenyl, benzyl, pyrrole, imidiz, pyridine, pyrimidine]

    global functional_groups_dict
    functional_groups_dict = {
        ester : 'ester',
        ether : 'ether',
        hydroxy : 'hydroxy',
        carbox_acid : 'carbox_acid',
        aldehyde :'aldehyde',
        anhydr : 'anhydr',
        amine : 'amine',
        amine_2 : 'amine_2',
        amine_3: 'amine_3',
        amide : 'amide',
        amide_2 : 'amide_2',
        amide_3: 'amide_3',
        nitro : 'nitro',
        imine: 'imine',
        f_hal :'f_hal',
        cl_hal : 'cl_hal',
        br_hal : 'br_hal',
        i_hal : 'i_hal',
        alkene : 'alkene',
        alkyne : 'alkyne',
        alkyne_term : 'alkyne_term',
        phenyl : 'phenyl',
        benzyl : 'benzyl',
        pyrrole : 'pyrrole',
        imidiz : 'imidiz',
        pyridine : 'pyridine',
        pyrimidine : 'pyrimidine'
    }

    global groups_to_numbers
    groups_to_numbers = {
        'ester' : 1,
        'ether' : 2,
        'hydroxy' : 3,
        'carbox_acid' : 4,
        'aldehyde' : 5,
        'anhydr' : 6,
        'amine' : 7,
        'amine_2' : 8,
        'amine_3': 9,
        'amide' : 10,
        'amide_2' : 11,
        'amide_3' : 12,
        'nitro' : 13,
        'imine' : 14,
        'f_hal' : 15,
        'cl_hal' : 16,
        'br_hal' : 17,
        'i_hal' : 18,
        'alkene' : 19,
        'alkyne' : 20,
        'alkyne_term' : 21,
        'phenyl' : 22,
        'benzyl' : 23,
        'pyrrole' : 24,
        'imidiz' : 25,
        'pyridine' : 26,
        'pyrimidine' : 27
    }

    global groups_dict
    groups_dict = {
        'ester': ['C(=O)-[NH]-C', 'C[C](=O)[O][C](=O)C'],
        'ether': ['C-S-C'],
        'hydroxy': ['C(=O)[OH]', 'C(=O)[H]', 'C(=O)-Cl'],
        'carbox_acid': ['C(=O)-[H]', 'C(=O)-Cl'],
        'aldehyde': ['C(=O)[OH]', 'C(=O)-Cl', 'C-[OH]', 'C-S-[H]', 'C(=[NH])[H]'],
        'anhydr': ['C(=O)-[NH]-C','C(=O)OC', 'C(=O)-[NH]-C'],
        'amine': ['CC', 'C-[N+](=O)-[O-]'],
        'amide': ['C(=O)[OH]', 'C(=O)[H]'],
        'amide_2': ['C(=O)OC'],
        'nitro': ['C-[NH2]'],
        'f_hal': ['Cl', 'Br', 'I'],
        'cl_hal': ['F', 'Br', 'I'],
        'br_hal': ['F', 'Cl', 'I'],
        'i_hal': ['F', 'Cl', 'Br'],
        'alkyne_term': ['C#N'], 
        'phenyl' : ['Cc1ccc([OH])cc1','Cc1ccc(OC)cc1', 'Cc1ccc(C)cc1', 'Cc1ccccc1'],
        'benzyl': ['CCc1ccc([OH])cc1','CCc1ccc(OC)cc1', 'CCc1ccc(C)cc1', 'CCc1ccccc1'],
        'pyrrole': ['[NH]1C=NC=C1'],
        'imidiz': ['[NH]1C=CC=C1'],
        'pyridine': ['n1cnccc1'],
        'pyrimidine': ['n1ccccc1']
    }
    return all_func_groups, type_dict, functional_groups, functional_groups_dict, groups_to_numbers, groups_dict

def group_idxes_from_mol(lig):
    """
    Create data path/dir and sub folders for each file type, return error if exists

    
    Parameters
    ----------
    lig : RDKIT molecule
        Ligand of interest.
    Returns
    -------
    match_indexes : dict
        Dictionary of functional groups and their corresponding atom indexes in the ligand.
    """
    match_indexes = {}
    mol = lig
    mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol))])))[1]
    mol_renum = Chem.RenumberAtoms(mol, mol_neworder)
    for j in functional_groups:
        k = Chem.MolFromSmarts(j)
        if mol_renum.HasSubstructMatch(k):
            idxes = mol_renum.GetSubstructMatches(k)
            idxes_list = list(idxes)
            for index in idxes:
                for subind in index:
                    if len(match_indexes[subind]) != 0:
                        temp_match = list(match_indexes[subind].values())
                        match_indexes[subind] = temp_match.append(str(functional_groups_dict[j]))
                    else:
                        match_indexes[subind] = list(str(functional_groups_dict[j]))
    return match_indexes
