import pandas as pd
import numpy as np
from config import *
from predict import ecif_predict
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


def my_test():
    predict_on_core('6.0')
    predict_on_core_ag(ecif_example)


if __name__ == '__main__':
    my_test()
