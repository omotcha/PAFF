import sklearn.metrics
from autogluon.tabular import TabularDataset, TabularPredictor
import os
from config import *
from util.ECIF import *
import csv
from DatasetMng import IndexMng
from autogluon.core.metrics import make_scorer


def get_ecif_dict():
    ret = {}
    with open('../Preprocess/BindingData.csv', 'r') as f:
        f.readline()
        reader = csv.reader(f)
        for row in reader:
            ret[row[0]] = row[2]
    return ret


def train_with_autogluon(distance_cutoffs, split_factor):
    """
    training tabular data (like ECIF_6.0.csv) using autogluon
    :param distance_cutoffs:
    :param split_factor: size(test) / size(test+train),
                         if split_factor not given, you should split it yourself
    :return:
    """
    for d in distance_cutoffs:
        print("\n")
        ecif = TabularDataset(os.path.join(tmpdata_dir, 'ecif_data', 'ECIF_{}.csv'.format(d)))
        ligand_descriptors = TabularDataset(os.path.join(tmpdata_dir, 'ecif_data', 'RDKit_Descriptors.csv'))
        binding_data = TabularDataset("../Preprocess/BindingData.csv")
        # Merge descriptors
        ecif = ecif.merge(ligand_descriptors, left_on="PDB", right_on="PDB")
        ecif = ecif.merge(binding_data, left_on="PDB", right_on="PDB")
        if not split_factor:
            # uses dataset split given by ECIF paper
            # Split training and test sets
            train_x = ecif[ecif["SET"] == "Train"][list(ecif.columns)[1:-2]]
            train_y = ecif[ecif["SET"] == "Train"]["pK"]

            test_x = ecif[ecif["SET"] == "Test"][list(ecif.columns)[1:-2]]
            test_y = ecif[ecif["SET"] == "Test"]["pK"]
        else:
            from sklearn.model_selection import train_test_split
            X = ecif[list(ecif.columns)[1:-2]]
            y = ecif['pK']
            train_x, test_x, train_y, test_y = train_test_split(
                X,
                y,
                test_size=split_factor if 0 < split_factor < 1 else 0.25,
                random_state=0)
        print("training...\n")
        ag_r2_scorer = make_scorer(name='r2',
                                   score_func=sklearn.metrics.r2_score,
                                   optimum=1,
                                   greater_is_better=True)
        predictor = TabularPredictor(label='pK',
                                     eval_metric='root_mean_squared_error').fit(train_data=train_x.join(train_y))
        print("\n - finished -\n")
        leaderboard = predictor.leaderboard(test_x.join(test_y), extra_metrics=[ag_r2_scorer], silent=True)
        print(leaderboard[['model', 'score_val', 'r2']])


def predict_with_autogluon_by_id(ids, d, model_dir, select_model=None):
    print("\n")
    predictor = TabularPredictor.load(model_dir)
    ecif_dict = get_ecif_dict()
    ecif_core_ids = IndexMng.get_index_from_dir(ecif_core)
    ecif_2016_refined_ids = IndexMng.get_index_from_dir(ecif_2016_refined)
    ecif_2019_refined_ids = IndexMng.get_index_from_dir(ecif_2019_refined)
    ecif_2019_general_minus_refined_ids = IndexMng.get_index_from_dir(ecif_2019_general_minus_refined)
    helper = ECIF(2016)
    pred_results = []
    actual_pk = []
    for iid in ids:
        protein_file = ""
        ligand_file = ""
        actual_pk.append(ecif_dict[iid])
        if iid in ecif_core_ids:
            protein_file = os.path.join(ecif_core, iid, "{}_protein.pdb".format(iid))
            ligand_file = os.path.join(ecif_core, iid, "{}_ligand.sdf".format(iid))

        elif iid in ecif_2016_refined_ids:
            protein_file = os.path.join(ecif_2016_refined, iid, "{}_protein.pdb".format(iid))
            ligand_file = os.path.join(ecif_2016_refined, iid, "{}_ligand.sdf".format(iid))

        elif iid in ecif_2019_general_minus_refined_ids:
            protein_file = os.path.join(ecif_2019_general_minus_refined, iid, "{}_protein.pdb".format(iid))
            ligand_file = os.path.join(ecif_2019_general_minus_refined, iid, "{}_ligand.sdf".format(iid))

        elif iid in ecif_2019_refined_ids:
            protein_file = os.path.join(ecif_2019_refined, iid, "{}_protein.pdb".format(iid))
            ligand_file = os.path.join(ecif_2019_refined, iid, "{}_ligand.sdf".format(iid))

        ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
        # elements_list = helper.get_elements(protein_file, ligand_file, float(d))
        ligand_descriptors = helper.get_ligand_features_by_file(ligand_file)
        ecif_list = pd.DataFrame([ecif_list], columns=helper.get_possible_ecif())
        ligand_descriptors = pd.DataFrame([ligand_descriptors], columns=helper.get_ligand_descriptors())
        Descriptors = ecif_list.join(ligand_descriptors)
        if select_model:
            pred = predictor.predict(Descriptors, model=select_model)
        else:
            pred = predictor.predict(Descriptors)
        pred_results.append(pred[0])

    prediction = pd.DataFrame(pred_results, columns=["Predicted pK"])
    actual = pd.DataFrame(actual_pk, columns=["Actual pK"])
    result = pd.DataFrame(ids, columns=["ID"]).join(prediction).join(actual)
    print(result)


def exp_step_1():
    train_with_autogluon(['6.0'], split_factor=False)


def exp_step_2():
    predict_with_autogluon_by_id(
        ids=['1a1b', '1a30'],
        d='6.0',
        model_dir=os.path.join(autogluon_model_dir, 'ag-20220629_071709'),
        select_model=False)


def my_test():
    exp_step_1()


if __name__ == '__main__':
    my_test()
