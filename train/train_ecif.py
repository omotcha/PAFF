import pandas as pd
import os
from config import *
import sklearn.tree
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from math import sqrt
import pickle


def train(distance_cutoffs, save_model=True):
    for d in distance_cutoffs:
        print("\n")
        ecif = pd.read_csv(os.path.join(tmpdata_dir, 'ecif_data', 'ECIF_{}.csv'.format(d)))
        ligand_descriptors = pd.read_csv(os.path.join(tmpdata_dir, 'ecif_data', 'RDKit_Descriptors.csv'))
        binding_data = pd.read_csv("../Preprocess/BindingData.csv")

        # Merge descriptors
        ecif = ecif.merge(ligand_descriptors, left_on="PDB", right_on="PDB")
        ecif = ecif.merge(binding_data, left_on="PDB", right_on="PDB")

        # Split training and test sets
        x_train = ecif[ecif["SET"] == "Train"][list(ecif.columns)[1:-2]]
        y_train = ecif[ecif["SET"] == "Train"]["pK"]

        x_test = ecif[ecif["SET"] == "Test"][list(ecif.columns)[1:-2]]
        y_test = ecif[ecif["SET"] == "Test"]["pK"]
        print("training...\n")
        RF = RandomForestRegressor(random_state=1206, n_estimators=500, n_jobs=8, oob_score=True, max_features=0.33)
        RF.fit(x_train, y_train)

        y_pred_RF = RF.predict(x_test)

        GBT = GradientBoostingRegressor(random_state=1206, n_estimators=20000, max_features="sqrt", max_depth=8,
                                        min_samples_split=3, learning_rate=0.005, loss="ls", subsample=0.7)
        GBT.fit(x_train, y_train)

        y_pred_GBT = GBT.predict(x_test)

        print("Pearson correlation coefficient for RF: ", pearsonr(y_test, y_pred_RF)[0])
        print("RMSE for RF:", sqrt(mean_squared_error(y_test, y_pred_RF)))
        print("Pearson correlation coefficient for GBT: ", pearsonr(y_test, y_pred_GBT)[0])
        print("RMSE for GBT:", sqrt(mean_squared_error(y_test, y_pred_GBT)))

        if save_model:
            print("Saving model...\n")
            pickle.dump(GBT, open(os.path.join(ecif_model_dir, 'ecif_gbt_{}.pkl'.format(d)), 'wb'))

        print("\n - finished -\n")


# def test():
#     train(['6.0'], False)


if __name__ == '__main__':
    train(['6.0'], save_model=True)
