"""Functions used to create and validate ligands from a SMILES string."""
from rdkit import Chem
from openbabel import pybel
import ipywidgets as widgets
from ipywidgets import Text, Layout, Label, Box, HBox, Output
from IPython.display import display

def create_ligands_from_smiles(num_of_ligs):
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
    form_item_layout = Layout(display='flex',flex_flow='row',justify_content='space-between')
    global names_for_ligs
    global smiles_for_ligs
    names_for_ligs = {}
    smiles_for_ligs = {}
    labels_name = {}
    labels_smiles = {}
    number = 1
    while number <= num_of_ligs.value:
        temp_name = "name" + str(number)
        temp_smile = "scratch" + str(number)
        names_for_ligs[temp_name] = Text(value = '', placeholder=f'Type the name of ligand {number} with no spaces', disabled=False)
        smiles_for_ligs[temp_smile] = Text(value = '', placeholder=f'Type in ligand {number} using SMILE codes', disabled=False)
        labels_name[temp_name] = Label(value=f'Name of ligand {number}')
        labels_smiles[temp_smile] = Label(value=f'SMILES string for ligand {number}')
        number += 1
    form_items1 = []
    form_items2 = []
    for list_number, name in enumerate(names_for_ligs):
        temp_box = "box" + str(list_number)
        new_number = list_number + 1
        form_items1.append(Box([labels_name["name" + str(new_number)], names_for_ligs[name]], layout=form_item_layout))
        
    for list_number, smiles in enumerate(smiles_for_ligs):
        temp_box = "box" + str(list_number)
        new_number = list_number + 1
        form_items2.append(Box([labels_smiles["scratch" + str(new_number)], smiles_for_ligs[smiles]], layout=form_item_layout))

    form1 = Box(form_items1, layout = Layout(
        display = 'flex',
        flex_flow = 'column',
        border = 'solid 2px',
        align_items = 'stretch',
        width = '50%'
    ))
    form2 = Box(form_items2, layout = Layout(
        display = 'flex',
        flex_flow = 'column',
        border = 'solid 2px',
        align_items = 'stretch',
        width = '50%'
    ))

    form = HBox([form1, form2])
    return form

def create_mols_from_smiles(num_of_ligs):
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
    name_vals = {}
    scratch_vals = {}
    for value in names_for_ligs.keys():
        name_vals[value] = names_for_ligs[value].value

    for value in smiles_for_ligs.keys():
        scratch_vals[value] = smiles_for_ligs[value].value
    
    smiles = []
    smile_names = []
    a = 0
    while a < num_of_ligs.value:
        name_temp = "name" + str(a + 1)
        scratch_temp = "scratch" + str(a + 1)
        lig_name = name_vals[name_temp]
        lig_scratch = scratch_vals[scratch_temp]
        lig_test = Chem.MolFromSmiles(lig_scratch)
        if (len(lig_scratch) < 2000) & (lig_test is not None):
            smile_names.append(lig_name)
            smiles.append(lig_scratch)
        a += 1   
    out=pybel.Outputfile(filename='data/MOL2_files/InputMols.mol2',format='mol2',overwrite=True)
    for index, smi in enumerate(smiles):
        mol = pybel.readstring(string=smi,format='smiles')
        mol.title= str(smile_names[index])
        mol.make3D('mmff94s')
        mol.localopt(forcefield = 'mmff94s', steps = 500)
        out.write(mol)
    out.close()
    return name_vals, scratch_vals

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    create_ligands_from_smiles()