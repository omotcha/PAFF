import pandas as pd
import numpy as np
from config import *
from DatasetMng.IndexMng import get_index_from_dir
import matplotlib.pyplot as plt
import pickle
from sklearn.metrics import r2_score, mean_squared_error

affinity_data = pd.read_csv('../Preprocess/BindingData.csv', comment='#')
ecif_data = pd.read_csv(os.path.join(ecif_data_dir, 'ECIF_6.0.csv'), comment='#')
ld_data = pd.read_csv(os.path.join(ecif_data_dir, 'RDKit_Descriptors.csv'), comment='#')


def query_pk_by_id(qid):
    return float(affinity_data[affinity_data['PDB'] == qid]['pK'])


def query_described_value(qid):
    ecif_cut = ecif_data[ecif_data['PDB'] == qid].iloc[:, 1:]
    ld_cut = ld_data[ld_data['PDB'] == qid].iloc[:, 1:]
    return ecif_cut.join(ld_cut)


def predict_on_core(d):
    core_ids = get_index_from_dir(ecif_core)
    actual = []
    predicted = []
    model = pickle.load(open(os.path.join(ecif_model_dir, 'ecif_gbt_{}.pkl'.format(d)), 'rb'))
    for id in core_ids:
        actual.append(query_pk_by_id(id))
        predicted.append(model.predict(query_described_value(id))[0])
    plt.title(label='ECIF::GBT')
    plt.scatter(actual, predicted, color="orange")
    # plt.plot([actual.min(), actual.max()], [predicted.min(), predicted.max()], 'r--', lw=2)
    plt.xlabel('actual')
    plt.ylabel('predicted')
    plt.xlim(0.0, 16.0)
    plt.ylim(0.0, 16.0)
    line_param = np.polyfit(actual, predicted, deg=1)
    y = [i*line_param[0]+line_param[1] for i in actual]
    plt.plot(actual, y, color="red")
    rmse = mean_squared_error(actual, predicted) ** 0.5
    r = r2_score(actual, predicted) ** 0.5
    plt.text(2, 14, 'rmse: {}'.format(rmse))
    plt.text(2, 12, 'r: {}'.format(r))

    plt.show()


def predict_on_core_ag(model):
    core_ids = get_index_from_dir(ecif_core)
    actual = []
    predicted = []
    model = pickle.load(open(model, 'rb'))
    for id in core_ids:
        actual.append(query_pk_by_id(id))
        predicted.append(model.predict(query_described_value(id))[0])
    plt.title(label='ECIF::Autogluon')
    plt.scatter(actual, predicted, color="orange")
    # plt.plot([actual.min(), actual.max()], [predicted.min(), predicted.max()], 'r--', lw=2)
    plt.xlabel('actual')
    plt.ylabel('predicted')
    plt.xlim(0.0, 16.0)
    plt.ylim(0.0, 16.0)
    line_param = np.polyfit(actual, predicted, deg=1)
    y = [i * line_param[0] + line_param[1] for i in actual]
    plt.plot(actual, y, color="red")
    rmse = mean_squared_error(actual, predicted) ** 0.5
    r = r2_score(actual, predicted) ** 0.5
    plt.text(2, 14, 'rmse: {}'.format(rmse))
    plt.text(2, 12, 'r: {}'.format(r))

    plt.show()


def predict_ex(f_pro, f_lig, d):
    """
    predict on protein-ligand pairs not in core set
    :param f_pro: protein file, .pdb
    :param f_lig: ligand file, .sdf
    :param d: using cutoff distance
    :return:
    """
    from util.ECIF import ECIF
    ecif_helper = ECIF(2016)
    ecif = ecif_helper.get_ecif(f_pro, f_lig, float(d))
    ld = ecif_helper.get_ligand_features_by_file(f_lig)
    model = pickle.load(open(os.path.join(ecif_model_dir, 'ecif_gbt_{}.pkl'.format(d)), 'rb'))
    data = ecif + list(ld)
    cols = ecif_helper.get_possible_ecif() + ecif_helper.get_ligand_descriptors()
    data_f = pd.DataFrame([data], columns=cols)
    return model.predict(data_f)[0]


def predict_ex_ag(f_pro, f_lig, model):
    """
    predict on protein-ligand pairs not in core set
    :param f_pro: protein file, .pdb
    :param f_lig: ligand file, .sdf
    :param model:
    :return:
    """
    from util.ECIF import ECIF
    ecif_helper = ECIF(2016)
    ecif = ecif_helper.get_ecif(f_pro, f_lig, float(6.0))
    ld = ecif_helper.get_ligand_features_by_file(f_lig)
    m = pickle.load(open(model, 'rb'))
    data = ecif + list(ld)
    cols = ecif_helper.get_possible_ecif() + ecif_helper.get_ligand_descriptors()
    data_f = pd.DataFrame([data], columns=cols)
    return m.predict(data_f)[0]


def my_test():
    predict_on_core('6.0')
    # predict_on_core_ag(ecif_example)


def estim_ext(ids, d):
    actual = []
    predicted = []
    model = pickle.load(open(os.path.join(ecif_model_dir, 'ecif_gbt_{}.pkl'.format(d)), 'rb'))
    for id in ids:
        actual.append(query_pk_by_id(id))
        predicted.append(model.predict(query_described_value(id))[0])
    print("\n")
    print("pred:")
    print(predicted)
    print("actual:")
    print(actual)


def estim_not_ext():
    print("\n")
    ids = ['6V1J', '6V1M', '6V1O']
    preds = []
    preds_ag = []
    for i in ids:
        f_prot = "C:\\Users\\user\\Desktop\\bak\\{}\\{}_protein.pdb".format(i, i)
        f_ligd = "C:\\Users\\user\\Desktop\\bak\\{}\\{}_ligand.sdf".format(i, i)
        preds.append(predict_ex(f_prot, f_ligd, '6.0'))
        preds_ag.append(predict_ex_ag(f_prot, f_ligd, ecif_example))

    print("\n")
    print(preds)
    print(preds_ag)
    for p in preds:
        print(10**(-p))


if __name__ == '__main__':
    print("\nnot existing")
    estim_not_ext()
    print("\nexisting")
    estim_ext(['2j4i', '2p4y'], '6.0')
    # my_test()
