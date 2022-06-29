import os

import pandas as pd
import pickle
from config import *
from DatasetMng import IndexMng
from util.ECIF import *
import csv


def get_ecif_dict():
    ret = {}
    with open('../Preprocess/BindingData.csv', 'r') as f:
        f.readline()
        reader = csv.reader(f)
        for row in reader:
            ret[row[0]] = row[2]
    return ret


def ecif_pred_by_id(ids, d):
    ecif_dict = get_ecif_dict()
    ecif_core_ids = IndexMng.get_index_from_dir(ecif_core)
    ecif_2016_refined_ids = IndexMng.get_index_from_dir(ecif_2016_refined)
    ecif_2019_refined_ids = IndexMng.get_index_from_dir(ecif_2019_refined)
    ecif_2019_general_minus_refined_ids = IndexMng.get_index_from_dir(ecif_2019_general_minus_refined)
    helper = ECIF(2016)
    pred_results = []
    actual_pk = []
    for id in ids:
        protein_file = ""
        ligand_file = ""
        actual_pk.append(ecif_dict[id])
        if id in ecif_core_ids:
            protein_file = os.path.join(ecif_core, id, "{}_protein.pdb".format(id))
            ligand_file = os.path.join(ecif_core, id, "{}_ligand.sdf".format(id))

        elif id in ecif_2016_refined_ids:
            protein_file = os.path.join(ecif_2016_refined, id, "{}_protein.pdb".format(id))
            ligand_file = os.path.join(ecif_2016_refined, id, "{}_ligand.sdf".format(id))

        elif id in ecif_2019_general_minus_refined_ids:
            protein_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_protein.pdb".format(id))
            ligand_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_ligand.sdf".format(id))

        elif id in ecif_2019_refined_ids:
            protein_file = os.path.join(ecif_2019_refined, id, "{}_protein.pdb".format(id))
            ligand_file = os.path.join(ecif_2019_refined, id, "{}_ligand.sdf".format(id))

        ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
        # elements_list = helper.get_elements(protein_file, ligand_file, float(d))
        ligand_descriptors = helper.get_ligand_features_by_file(ligand_file)
        ecif_list = pd.DataFrame([ecif_list], columns=helper.get_possible_ecif())
        ligand_descriptors = pd.DataFrame([ligand_descriptors], columns=helper.get_ligand_descriptors())
        Descriptors = ecif_list.join(ligand_descriptors)
        model = pickle.load(open(os.path.join(ecif_model_dir, 'ecif_gbt_{}.pkl'.format(d)), 'rb'))
        pred_results.append(model.predict(Descriptors)[0])

    prediction = pd.DataFrame(pred_results, columns=["Predicted pK"])
    actual = pd.DataFrame(actual_pk, columns=["Actual pK"])
    result = pd.DataFrame(ids, columns=["ID"]).join(prediction).join(actual)
    print(result)


if __name__ == '__main__':
    ecif_pred_by_id(['1a1b', '1a30'], '6.0')
