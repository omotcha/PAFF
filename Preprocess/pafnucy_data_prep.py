import numpy as np
import pandas as pd
import h5py
import os
import re
import csv
from tqdm import tqdm
import matplotlib.pyplot as plt
from config import *
from openbabel import pybel
from tfbioSupport.tfbio.data import Featurizer

import warnings

import seaborn as sns

# omotcha: base dir of dataset on Windows
dataset_dir = 'D:\\AlexYu\\datasets\\dataset\\pdbbind2016'
dataset2013_dir = 'D:\\AlexYu\\datasets\\dataset\\pdbbind2013'


def get_aff_data():
    general_index = os.path.join(dataset_dir,
                                 'PDBbind_2016_plain_text_index', 'index', 'INDEX_general_PL_data.2016')

    aff_data = []
    with open(general_index, "r") as f:
        lines = f.readlines()
        for line in lines:
            seg = re.split(' +', line)
            if seg[0] != '#':
                aff_data.append([seg[0], seg[3]])
    return aff_data


def write_aff_data():
    with open(os.path.join(tmpdata_dir, 'affinity_data_1.csv'), "w") as f:
        f.write('pdbid,-logKd/Ki\n')
        for aff in get_aff_data():
            f.write('{},{}\n'.format(aff[0], aff[1]))


def get_missing():
    aff_data = get_aff_data()
    # omotcha: this aff_data should be initialized first by get_affinity_data
    exclusion_dir = os.path.join(dataset_dir,
                                 'general-set-except-refined')
    refined_dir = os.path.join(dataset_dir,
                               'refined-set')
    exclusion_ids = os.listdir(exclusion_dir)
    refined_ids = os.listdir(refined_dir)
    missing = [d[0] for d in aff_data if d[0] not in exclusion_ids and d[0] not in refined_ids]
    return missing


def get_affinity_data():
    missing = get_missing()
    affinity_data = pd.read_csv(os.path.join(tmpdata_dir, 'affinity_data_1.csv'), comment='#')
    affinity_data = affinity_data[~np.in1d(affinity_data['pdbid'], list(missing))]
    # affinity data should not contain NaN
    assert not affinity_data['-logKd/Ki'].isnull().any()
    return affinity_data


def get_cleaned_affinity_data():
    # cleaned datafile can be created by calling set_separation_and_exclude_2013()
    assert os.path.exists(os.path.join(tmpdata_dir, 'affinity_data_cleaned_1.csv'))
    affinity_data = pd.read_csv(os.path.join(tmpdata_dir, 'affinity_data_cleaned_1.csv'), comment='#')
    # affinity data should not contain NaN
    assert not affinity_data['-logKd/Ki'].isnull().any()
    return affinity_data


def set_separation_and_exclude_2013():
    affinity_data = get_affinity_data()

    # omotcha: note that refined set is trimmed cuz some proteins cannot be protonized and charged later,
    #          so it cannot be judged from the index file, it should be judged directly from the refined set,
    #          it is the same with core set
    # refined_index = os.path.join(dataset_dir,
    #                              'PDBbind_2016_plain_text_index', 'index', 'INDEX_refined_data.2016')
    # refined_set = set()
    # with open(refined_index, "r") as f:
    #     lines = f.readlines()
    #     for line in lines:
    #         seg = re.split(' +', line)
    #         if seg[0] != '#':
    #             refined_set.add(seg[0])
    refined_dir = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\refined-set\\"
    refined_list = [k for k in os.listdir(refined_dir)
                    if len(k) == 4 and os.path.isfile(refined_dir + "%s\\%s_pocket.mol2" % (k, k))]
    refined_set = set(refined_list)

    core_dir = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\coreset\\"
    core_list = [k for k in os.listdir(core_dir)
                 if len(k) == 4 and os.path.isfile(core_dir + "%s\\%s_pocket.mol2" % (k, k))]
    core_set = set(core_list)

    general_set = set(affinity_data['pdbid'])

    # core_index = os.path.join(dataset_dir,
    #                           'PDBbind_2016_plain_text_index', 'index', 'INDEX_core_data.2016')
    # core_set = set()
    # if not os.path.exists(core_index):
    #     print('Creating index file of core set:\n')
    #     core_dir = os.path.join(dataset_dir, 'coreset')
    #     core_ids = os.listdir(core_dir)
    #     core_index_data = []
    #     with open(refined_index, "r") as f:
    #         reader = csv.reader(f)
    #         for row in reader:
    #             if row[0][0:4] in core_ids:
    #                 core_index_data.append(row)
    #     with open(core_index, "w") as f:
    #         for row in core_index_data:
    #             f.writelines(str(row[0])+'\n')
    #     print('Done.\n')
    #
    # with open(core_index, "r") as f:
    #     lines = f.readlines()
    #     for line in lines:
    #         seg = re.split(' +', line)
    #         if seg[0] != '#':
    #             core_set.add(seg[0])

    # core set should be only included from refined set
    assert core_set & refined_set == core_set
    # refined set should be only included from general set
    assert refined_set & general_set == refined_set
    core2013_dir = os.path.join(dataset2013_dir, 'coreset')
    core2013_ids = os.listdir(core2013_dir)
    core2013 = set(core2013_ids)
    affinity_data['include'] = True
    affinity_data.loc[np.in1d(affinity_data['pdbid'], list(core2013 & (general_set - core_set))), 'include'] = False
    affinity_data.loc[np.in1d(affinity_data['pdbid'], list(general_set)), 'set'] = 'general'
    affinity_data.loc[np.in1d(affinity_data['pdbid'], list(refined_set)), 'set'] = 'refined'
    affinity_data.loc[np.in1d(affinity_data['pdbid'], list(core_set)), 'set'] = 'core'
    print(affinity_data[affinity_data['include']].groupby('set').apply(len).loc[['general', 'refined', 'core']])
    affinity_data[['pdbid']].to_csv('pdb_1.ids', header=False, index=False)
    affinity_data[['pdbid', '-logKd/Ki', 'set', 'include']].\
        to_csv(os.path.join(tmpdata_dir, 'affinity_data_cleaned_1.csv'), index=False)
    return general_set, refined_set, core_set, core2013


def prepare_pockets(file_path):
    pdb_list = [k for k in os.listdir(file_path) if
                len(k) == 4 and not os.path.isfile(file_path + "%s/%s_pocket.mol2" % (k, k))]
    print(len(pdb_list))
    count = 0
    for name in tqdm(pdb_list):
        if len(name) != 4:
            continue
        PDB_file = file_path + name + '/' + name
        os.system(" chimera --nogui process.py %s" % PDB_file)
        count += 1
    print(count)


def prepare_pockets_all():
    pass


def parse_molecule():
    affinity_data = get_cleaned_affinity_data()
    # omotcha: coreset used for test, but why here name 'refined-set' as 'core'?
    dataset_path = {'general': 'general-set-except-refined', 'refined': 'refined-set', 'core': 'coreset'}
    featurizer = Featurizer()
    charge_idx = featurizer.FEATURE_NAMES.index('partialcharge')

    with h5py.File('%s\\core2013.hdf' % dataset_dir, 'w') as g:
        j = 0

        for dataset_name, data in affinity_data.groupby('set'):

            print(dataset_name, 'set')
            i = 0
            ds_path = dataset_path[dataset_name]

            with h5py.File('%s\\%s.hdf' % (dataset_dir, dataset_name), 'w') as f:
                print("creating {}.hdf".format(dataset_name))
                for _, row in data.iterrows():

                    name = row['pdbid']
                    affinity = row['-logKd/Ki']

                    ligand = next(pybel.readfile('mol2', '%s\\%s\\%s\\%s_ligand.mol2' % (dataset_dir, ds_path, name, name)))
                    # do not add the hydrogens! they are in the strucutre and it would reset the charges

                    try:
                        pocket = next(pybel.readfile('mol2', '%s\\%s\\%s\\%s_pocket.mol2' % (dataset_dir, ds_path, name, name)))
                        # do not add the hydrogens! they were already added in chimera and it would reset the charges
                    except:
                        warnings.warn('no pocket for %s (%s set)' % (name, dataset_name))
                        continue

                    ligand_coords, ligand_features = featurizer.get_features(ligand, molcode=1)
                    assert (ligand_features[:, charge_idx] != 0).any()
                    pocket_coords, pocket_features = featurizer.get_features(pocket, molcode=-1)
                    assert (pocket_features[:, charge_idx] != 0).any()

                    centroid = ligand_coords.mean(axis=0)
                    ligand_coords -= centroid
                    pocket_coords -= centroid

                    data = np.concatenate((np.concatenate((ligand_coords, pocket_coords)),
                                           np.concatenate((ligand_features, pocket_features))), axis=1)

                    if row['include']:
                        dataset = f.create_dataset(name, data=data, shape=data.shape, dtype='float32',
                                                   compression='lzf')
                        dataset.attrs['affinity'] = affinity
                        i += 1
                    else:
                        dataset = g.create_dataset(name, data=data, shape=data.shape, dtype='float32',
                                                   compression='lzf')
                        dataset.attrs['affinity'] = affinity
                        j += 1

            print('prepared', i, 'complexes')
        print('excluded', j, 'complexes')


def write_protein_data():
    protein_data = pd.read_csv('%s\\PDBbind_2016_plain_text_index\\index\\INDEX_general_PL_name.2016' % dataset_dir,
                               comment='#', sep='  ', engine='python', na_values='------',
                               header=None, names=['pdbid', 'year', 'uniprotid', 'name'])

    # print(protein_data.head())
    # we assume that PDB IDs are unique
    assert ~protein_data['pdbid'].duplicated().any()
    protein_data = protein_data[np.in1d(protein_data['pdbid'], get_cleaned_affinity_data()['pdbid'])]

    # check for missing values
    # print(protein_data.isnull().any())
    # print(protein_data[protein_data['name'].isnull()])

    # fix rows with wrong separators between protein ID and name

    for idx, row in protein_data[protein_data['name'].isnull()].iterrows():
        uniprotid = row['uniprotid'][:6]
        name = row['uniprotid'][7:]
        protein_data.loc[idx, ['uniprotid', 'name']] = [uniprotid, name]

    # print(protein_data.isnull().any())

    protein_data.to_csv(os.path.join(tmpdata_dir, 'protein_data_1.csv'), index=False)


def exp_steps_1():
    write_aff_data()


def exp_steps_2():
    print(len(get_missing()))


def exp_steps_3():
    set_separation_and_exclude_2013()


def exp_steps_4():
    parse_molecule()


def exp_steps_5():
    write_protein_data()


def test():
    print('\n')
    exp_steps_5()


if __name__ == '__main__':
    test()