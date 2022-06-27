import os

import numpy as np
import pandas as pd
import h5py
import tfbioSupport
from config import *


def analysis(path):
    charge_std = 0.41113196291076504
    featurizer = tfbioSupport.tfbio.data.Featurizer()
    max_dist = 10
    box_size = 21
    columns = {name: i for i, name in enumerate(featurizer.FEATURE_NAMES)}
    num_features = len(columns)
    affinity = []
    coords = []
    features = []
    ids = []

    with h5py.File(os.path.join(base_dir, 'core2013.hdf'), 'r') as f:
        for pdb_id in f:
            ids.append(pdb_id)
            dataset = f[pdb_id]
            coords.append(dataset[:, :3])
            features.append(dataset[:, 3:])

            affinity.append(dataset.attrs['affinity'])

    affinity = np.reshape(affinity, (-1, 1))
    # print(affinity)
    batch_grid = []

    for crd, f in zip(coords, features):
        batch_grid.append(tfbioSupport.tfbio.data.make_grid(crd, f))

    batch_grid = np.vstack(batch_grid)
    batch_grid[..., columns['partialcharge']] /= charge_std


if __name__ == '__main__':
    # analysis("D:\\AlexYu\\datasets\\dataset\\2016_pred\\output-2022-06-23T13_02_58-predictions.csv")
    pass
