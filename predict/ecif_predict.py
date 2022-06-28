import os

import pandas as pd

from config import *
from DatasetMng import IndexMng
from util.ECIF import *
import csv

def get_ecif_ids():
    ids = []
    with open('../Preprocess/BindingData.csv', 'r') as f:
        f.readline()
        reader = csv.reader(f)
        for row in reader:
            ids.append(row[0])
    return ids


def ecif_pred_by_id(id, d):
    ecif_ids = get_ecif_ids()
    ecif_core_ids = IndexMng.get_index_from_dir(ecif_core)
    ecif_2016_refined_ids = IndexMng.get_index_from_dir(ecif_2016_refined)
    ecif_2019_refined_ids = IndexMng.get_index_from_dir(ecif_2019_refined)
    ecif_2019_general_minus_refined_ids = IndexMng.get_index_from_dir(ecif_2019_general_minus_refined)
    helper = ECIF(2016)
    if id in ecif_core_ids:
        protein_file = os.path.join(ecif_core, id, "{}_protein.pdb".format(id))
        ligand_file = os.path.join(ecif_core, id, "{}_ligand.sdf".format(id))
        ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
        # elements_list = helper.get_elements(protein_file, ligand_file, float(d))
        ligand_descriptors = helper.get_ligand_features_by_file(ligand_file)


    elif id in ecif_2016_refined_ids:
        protein_file = os.path.join(ecif_2016_refined, id, "{}_protein.pdb".format(id))
        ligand_file = os.path.join(ecif_2016_refined, id, "{}_ligand.sdf".format(id))
        ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
        # elements_list = helper.get_elements(protein_file, ligand_file, float(d))
        ligand_descriptors = helper.get_ligand_features_by_file(ligand_file)


    elif id in ecif_2019_general_minus_refined_ids:
        protein_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_protein.pdb".format(id))
        ligand_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_ligand.sdf".format(id))
        ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
        # elements_list = helper.get_elements(protein_file, ligand_file, float(d))
        ligand_descriptors = helper.get_ligand_features_by_file(ligand_file)
        ecif_list = pd.DataFrame(ecif_list, columns=helper.get_possible_ecif())
        ligand_descriptors = pd.DataFrame(ligand_descriptors, columns=helper.get_ligand_descriptors())

    elif id in ecif_2019_refined_ids:
        protein_file = os.path.join(ecif_2019_refined, id, "{}_protein.pdb".format(id))
        ligand_file = os.path.join(ecif_2019_refined, id, "{}_ligand.sdf".format(id))
        ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
        # elements_list = helper.get_elements(protein_file, ligand_file, float(d))
        ligand_descriptors = helper.get_ligand_features_by_file(ligand_file)



if __name__ == '__main__':
    ecif_pred_by_id('1a1b', '6.0')
