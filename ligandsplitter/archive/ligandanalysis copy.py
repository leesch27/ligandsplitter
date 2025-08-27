"""Provide the primary functions."""
import sys, os
from rdkit import Chem
import pandas as pd
from sklearn.model_selection import RandomizedSearchCV, train_test_split, cross_val_score, cross_validate
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.compose import ColumnTransformer, make_column_transformer
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

def get_vars():
    """
    Initialize essential variables for ligand analysis

    Parameters
    ----------
    None

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
