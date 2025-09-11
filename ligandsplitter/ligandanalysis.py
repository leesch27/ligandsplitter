"""Provide the primary functions."""
import sys, os
import numpy as np
import re
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Draw, Crippen, Lipinski
import pandas as pd
from sklearn.model_selection import RandomizedSearchCV, train_test_split, cross_val_score, cross_validate
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.compose import ColumnTransformer, make_column_transformer
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from .basefunctions import LigandVariables

vars = LigandVariables()

def group_idxes_from_mol(lig, renumber = False):
    """
    Get atom indices of functional groups in a ligand molecule.

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
    if renumber:
        mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol))])))[1]
        mol_renum = Chem.RenumberAtoms(mol, mol_neworder)
    else:
        mol_renum = mol
    for j in vars.functional_groups:
        k = Chem.MolFromSmarts(j)
        if mol_renum.HasSubstructMatch(k):
            idxes = mol_renum.GetSubstructMatches(k)
            dict_keys = list(match_indexes.keys())
            for index in idxes:
                for subind in index:
                    subind_string = str(subind)
                    if subind_string in dict_keys:
                        temp_match = match_indexes[subind_string]
                        match_indexes[subind_string] = temp_match.append(vars.functional_groups_dict[j])
                    else:
                        match_indexes[subind_string] = [vars.functional_groups_dict[j]]
    return match_indexes

def rf_classifier(data, method = ""):
    """
    Determine the importance of ligand features in whether they are orally bioactive or not using a Random Forest Classifier.
    
    Parameters
    ----------
    data : dataframe
        Data containing physical properties of ligand of interest.
    method : String
        Method to use for analysis. Options are "LRO5", "Ghose", or "Veber".

    Returns
    -------
    imps : dataframe
        Dataframe containing feature importances of the model.
    """

    # get feature and target variables
    if method == "LRO5":
        features = data.drop(columns = ["filename_hydrogens", "smiles", "mol_refractivity", "rotatable_bonds", "polar_surface_area", "orally_bioactive", "mol"])
    elif method == "Ghose":
        features = data.drop(columns = ["filename_hydrogens", "smiles", "rotatable_bonds", "polar_surface_area", "orally_bioactive", "mol"])
    elif method == "Veber":
        features = data.drop(columns = ["filename_hydrogens", "smiles", "orally_bioactive", "mol"])
    else:
        print("Please provide a valid method: LRO5, Ghose, or Veber.")
    target = data["orally_bioactive"]

    # define training and testing data
    X_train, X_test, y_train, y_test = train_test_split(features, target, test_size = 0.2)
    rf_c = RandomForestClassifier()
    rf_c.fit(X_train, y_train)
    scores_rf_c = cross_validate(rf_c, X_train, y_train, return_train_score=True)
    print(f"Initial cross-validation fit time: {scores_rf_c['fit_time']}")
    print(f"Initial cross-validation score time: {scores_rf_c['score_time']}")
    print(f"Initial cross-validation training scores: {scores_rf_c['train_score']}")
    print(f"Initial cross-validation testing scores: {scores_rf_c['test_score']}")

    # hyperparameter optimization
    print("Starting hyperparameter optimization...")
    rf_param_grid = {
        "max_depth": [1, 5, 10, 15, 20],
        "max_features": [1, 5, 10, 15, 20],
        "min_samples_split": [10, 20, 30, 40, 50],
        "min_samples_leaf": [5, 10, 15, 20]
    }
    rf_random_search = RandomizedSearchCV(
        RandomForestClassifier(), param_distributions=rf_param_grid, n_jobs=-1, n_iter=10, cv=5, random_state=123
    )
    print("Done!")

    # create and deploy optimized model
    print("Fitting optimized model...")
    rf_random_search.fit(X_train, y_train)
    optimized_rf_c = pd.DataFrame(rf_random_search.cv_results_)[["mean_test_score","param_max_depth","param_max_features", "param_min_samples_split", "param_min_samples_leaf", "mean_fit_time","rank_test_score",]].set_index("rank_test_score").sort_index().T
    max_depth_val = float(optimized_rf_c.iloc[1,1])
    max_features_val = float(optimized_rf_c.iloc[1,2])
    min_sample_split_val = float(optimized_rf_c.iloc[1,3])
    min_samples_leaf_val = float(optimized_rf_c.iloc[1,4])
    rf_optimal = RandomForestClassifier(max_depth = max_depth_val, max_features = max_features_val, min_samples_split = min_sample_split_val, min_samples_leaf = min_samples_leaf_val)
    rf_optimal.fit(X_train, y_train)
    rf_optimal.score(X_test, y_test)

    # obtain feature importance
    feature_names = list(X_train.columns)
    data = {
        "Importance": rf_optimal.feature_importances_,
    }
    imps = pd.DataFrame(data=data, index=feature_names,).sort_values(by="Importance", ascending=False)[:10]
    return imps

def rf_regressor(data):
    """
    Determine the importance of ligand features in determining binding affinity.
    
    Parameters
    ----------
    data : dataframe
        Data containing docking information between ligand/s of interest and receptor.

    Returns
    -------
    imps : dataframe
        Dataframe containing feature importances of the model.
    """

    features = data.drop(columns = ["Score"])
    target = data["Score"]
    X_train, X_test, y_train, y_test = train_test_split(features, target, test_size = 0.2)
    
    # initial training
    rf_r = RandomForestClassifier()
    rf_r.fit(X_train, y_train)
    scores_rf_r = cross_validate(rf_r, X_train, y_train, return_train_score=True)
    print(f"Initial cross-validation fit time: {scores_rf_r['fit_time']}")
    print(f"Initial cross-validation score time: {scores_rf_r['score_time']}")
    print(f"Initial cross-validation training scores: {scores_rf_r['train_score']}")
    print(f"Initial cross-validation testing scores: {scores_rf_r['test_score']}")

    # hyperparameter optimization
    print("Starting hyperparameter optimization...")
    rf_param_grid = {
        "max_depth": [1, 5, 10, 15, 20],
        "max_features": [1, 5, 10, 15, 20],
        "min_samples_split": [10, 20, 30, 40, 50],
        "min_samples_leaf": [5, 10, 15, 20]
    }
    rf_random_search = RandomizedSearchCV(RandomForestRegressor(), param_distributions=rf_param_grid, n_jobs=-1, n_iter=10, cv=5, random_state=123)
    print("Done!")

    # create and deploy optimized model
    print("Fitting optimized model...")
    rf_random_search.fit(X_train, y_train)
    optimized_rf_r = pd.DataFrame(rf_random_search.cv_results_)[["mean_test_score","param_max_depth","param_max_features", "param_min_samples_split", "param_min_samples_leaf", "mean_fit_time","rank_test_score",]].set_index("rank_test_score").sort_index().T
    max_depth_val = float(optimized_rf_r.iloc[1,1])
    max_features_val = float(optimized_rf_r.iloc[1,2])
    min_sample_split_val = float(optimized_rf_r.iloc[1,3])
    min_samples_leaf_val = float(optimized_rf_r.iloc[1,4])
    rf_optimal = RandomForestClassifier(max_depth = max_depth_val, max_features = max_features_val, min_samples_split = min_sample_split_val, min_samples_leaf = min_samples_leaf_val)
    rf_optimal.fit(X_train, y_train)
    rf_optimal.score(X_test, y_test)

    # obtain feature importance
    feature_names = list(X_train.columns)
    data = {
        "Importance": rf_optimal.feature_importances_,
    }
    imps = pd.DataFrame(data=data, index=feature_names,).sort_values(by="Importance", ascending=False)[:10]
    return imps

def number_of_atoms(atom_list, df):
    """
    Determine the number of atoms each ligand has for a given dataframe.
    
    Parameters
    ----------
    atom_list : list
        List of atom name abbreviations.
    df : dataframe
        Data containing initial data on ligand properties.

    Returns
    -------
    None
    """
    # determine the number of different heavy atoms in each ligand
    for i in atom_list:
        substruct_list = []
        for index, row in df.iterrows():
            smile_string = row['smiles']
            if len(i) == 1:
                string_finder_lower = re.findall(r'{}(?![aelu+][+\d])(?!([aeolu]+[+\d]))'.format(i.lower()), smile_string)
                string_finder_upper = re.findall(r'{}(?![aelu+][+\d])(?!([aeolu]+[+\d]))'.format(i), smile_string)
                substruct_list.append(len(string_finder_lower) + len(string_finder_upper))
            else:
                string_finder_brackets = re.findall(r'[\[]{}[\]]'.format(i), smile_string)
                string_finder_charged = re.findall(r'[\[]{}[+][+\d]'.format(i), smile_string)
                substruct_list.append(len(string_finder_brackets) + len(string_finder_charged))
        df['num_of_{}_atoms'.format(i)] = substruct_list

def atom_weights(df):
    """
    Calculate the weight of each ligand in a dataframe.
    
    Parameters
    ----------
    df : dataframe
        Dataframe of ligand properties that includes atom counts (number_of_atoms function).

    Returns
    -------
    None
    """
    # calculate weight of ligands
    global atom_weights_dict
    atom_weights_dict = {
        'C':12.0096,
        'N': 14.006,
        'O': 15.999,
        'F': 18.998,
        'Al': 26.981,
        'P': 30.974,
        'S': 32.059,
        'Cl': 35.45,
        'Cr': 51.9961,
        'Mn': 54.938,
        'Fe': 55.845,
        'Co': 58.933,
        'Ni': 58.693,
        'Cu': 63.546,
        'Zn': 65.38,
        'Ga': 69.723,
        'Ge': 72.630,
        'As': 74.921,
        'Br': 79.901,
        'Zr': 91.224,
        'Mo': 95.95,
        'Pd': 106.42,
        'Ag': 107.8682,
        'Cd': 112.414,
        'In': 114.818,
        'Sn': 118.71,
        'Sb': 121.760,
        'I': 126.904,
        'Ir': 192.217,
        'Pt': 195.08,
        'Au': 196.966570,
        'Hg': 200.592,
        'Pb': 207.2,
        'Bi': 208.980
    }
    ligand_weights = []
    for index, row in df.iterrows():
        ligand_atom_nums = sum(row[5:])
        weight_da = 0
        if row['num_of_heavy_atoms'] == ligand_atom_nums:
            for num, column in enumerate(row[5:]):
                column_title = list(df)[num + 5]
                atom_name = re.split("_", column_title)
                atom_type_weight = atom_weights_dict[atom_name[2]]
                weight_da = weight_da + (atom_type_weight *  column)
        weight_da = weight_da + ((row.iloc[3] - row.iloc[4]) * 1.007)
        ligand_weights.append(weight_da)
    df.insert(2, "molecular_weight", ligand_weights)

def chemical_physical_properties(df):
    """
    Calculate various properties for each ligand in a dataframe.
    
    Parameters
    ----------
    df : dataframe
        Data containing initial data on ligand properties.

    Returns
    -------
    None
    """
    # calculate logP (partition coefficient), hydrogen bond donors, hydrogen bond acceptors,
    # molar refractivity (Ghose filter), number of rotatable bonds (veber's rule) and polar surface
    # area (veber's rule) of ligands
    log_P = []
    H_donors = []
    H_acceptors = []
    mol_mr = []
    mol_rotatable = []
    tpsas = []
    for index, row in df.iterrows():
        mol = row.iloc[3]
        if type(mol) != float:
            log = Crippen.MolLogP(mol)
            log_P.append(log)
            donor = Lipinski.NumHDonors(mol)
            H_donors.append(donor)
            acceptor = Lipinski.NumHAcceptors(mol)
            H_acceptors.append(acceptor)
            mr = Crippen.MolMR(mol)
            mol_mr.append(mr)
            rotatable = Lipinski.NumRotatableBonds(mol)
            mol_rotatable.append(rotatable)
            psa = Descriptors.TPSA(mol)
            tpsas.append(psa)
        else:
            pass
    df.insert(3, "log_P", log_P)
    df.insert(4, "H_donors", H_donors)
    df.insert(5, "H_acceptors", H_acceptors)
    df.insert(6, "mol_refractivity", mol_mr)
    df.insert(7, "rotatable_bonds", mol_rotatable)
    df.insert(8, "polar_surface_area", tpsas)

def get_ligand_properties(lig_df):
    """
    Determine the importance of ligand features in determining binding affinity.
    
    Parameters
    ----------
    lig_df : dataframe
        Data containing initial data on ligand properties.

    Returns
    -------
    updated_df : dataframe
        Dataframe containing calculated ligand properties of interest.
    """
    # determine and record the number of atoms and the number of heavy atoms in each ligand
    global atom_abbriv
    atom_abbriv = ['C','N','O','F','Al','P','S','Cl','Cr','Mn','Fe','Co','Ni','Cu',
               'Zn','Ga','Ge','As','Br','Zr','Mo','Pd','Ag','Cd','In','Sn','Sb',
               'I','Ir','Pt','Au','Hg','Pb','Bi']

    mol_format = []
    atom_total = []
    atom_total_heavy = []
    updated_df = lig_df.copy()
    for index, row in updated_df.iterrows():
        try:
            mol = Chem.MolFromMol2File(row['filename_hydrogens'],sanitize=False)
            if mol is not None:
                mol_H = Chem.AddHs(mol)
                mol_format.append(mol_H)
                mol_atoms = mol_H.GetNumAtoms()
                atom_total.append(mol_atoms)
                mol_atoms_heavy = mol_H.GetNumHeavyAtoms()
                atom_total_heavy.append(mol_atoms_heavy)
            else:
                #currently only works for molecules containing only atoms with single letter names, need to fix
                string = row['smiles']
                string_alpha = re.findall(r'[a-zA-Z]', string)
                string_H = re.findall(r'[H]', string)
                mol_format.append(np.nan)
                atom_total.append(len(string_alpha))
                atom_total_heavy.append(len(string_alpha) - len(string_H))
        except OSError:
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is not None:
                mol_H = Chem.AddHs(mol)
                mol_format.append(mol_H)
                mol_atoms = mol_H.GetNumAtoms()
                atom_total.append(mol_atoms)
                mol_atoms_heavy = mol_H.GetNumHeavyAtoms()
                atom_total_heavy.append(mol_atoms_heavy)
            else:
                #currently only works for molecules containing only atoms with single letter names, need to fix
                string = row['smiles']
                string_alpha = re.findall(r'[a-zA-Z]', string)
                string_H = re.findall(r'[H]', string)
                mol_format.append(np.nan)
                atom_total.append(len(string_alpha))
                atom_total_heavy.append(len(string_alpha) - len(string_H))
    updated_df['mol'] = mol_format
    updated_df['num_of_atoms'] = atom_total
    updated_df['num_of_heavy_atoms'] = atom_total_heavy
    number_of_atoms(atom_abbriv, updated_df)
    atom_weights(updated_df)
    chemical_physical_properties(updated_df)
    
    return updated_df