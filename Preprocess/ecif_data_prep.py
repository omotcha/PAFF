import pandas as pd

from config import *
import os
import csv
from DatasetMng import IndexMng
from util.ECIF import *
from tqdm import tqdm

# make sure affinity data exists, it can be created by util/affinityDataMng
affinity_data = pd.read_csv(os.path.join(tmpdata_dir, 'affinity_data_1.csv'), comment='#')


def get_ecif_dict():
    ret = {}
    with open('BindingData.csv', 'r') as f:
        f.readline()
        reader = csv.reader(f)
        for row in reader:
            ret[row[0]] = row[2]
    return ret


def get_all_ign_index():
    """
    ign index folder has 3 index text files, named (valid|train|test).index.txt, containing all the protein indices
    :return:
    """
    with open(os.path.join(ign_index_dir, 'train.index.txt'), 'r') as ftr:
        tr = [i[0:4] for i in ftr.readlines()]
    with open(os.path.join(ign_index_dir, 'test.index.txt'), 'r') as fte:
        te = [i[0:4] for i in fte.readlines()]
    return tr + te


def get_all_ign_index_2013(test_num):
    def get_2013_core_ids():
        ids = []
        files = os.listdir(core2013_dir)
        for f in files:
            ids.append(f)
        return ids
    ids_2013 = get_2013_core_ids()
    with open(os.path.join(ign_index_dir, 'train.index.txt'), 'r') as ftr:
        tr = [i[0:4] for i in ftr.readlines()]
    count = 0
    te = []
    with open(os.path.join(ign_index_dir, 'test.index.txt'), 'r') as fte:
        for i in fte.readlines():
            if i[0:4] in ids_2013 and count < test_num:
                te.append(i[0:4])
                count = count + 1
                continue
            tr.append(i[0:4])
    return tr + te


def collect_ecif(distance_cutoffs, indices=None, tag=None):
    """
    collect ECIF and ELEMENTS data from pdb files and sdf files
    :param distance_cutoffs:
    :param indices: given ids for extra experiments,
                    note that currently this set of index should be subset of that of ECIF,
                    otherwise the not-included part would not be collected
    :param tag:     for extra experiments, to make a distinction between files created by ECIF
    :return:
    """
    if indices is not None:
        ecif_ids = indices
    else:
        ecif_dict = get_ecif_dict()
        ecif_ids = ecif_dict.keys()
    ecif_core_ids = IndexMng.get_index_from_dir(ecif_core)
    ecif_2016_refined_ids = IndexMng.get_index_from_dir(ecif_2016_refined)
    ecif_2019_refined_ids = IndexMng.get_index_from_dir(ecif_2019_refined)
    ecif_2019_general_minus_refined_ids = IndexMng.get_index_from_dir(ecif_2019_general_minus_refined)
    helper = ECIF(2016)
    print("\nCollecting ECIF and ELEMENTS data: ")
    not_listed_ids = set()
    for d in distance_cutoffs:
        print("\n distance_cutoff: {}\n".format(d))
        if tag is not None:
            f_ecif = open(os.path.join(tmpdata_dir, 'ecif_data', 'ECIF_{}_{}.csv'.format(d, tag)), 'a', newline='')
            f_elements = open(os.path.join(tmpdata_dir, 'ecif_data', 'ELEMENTS_{}_{}.csv'.format(d, tag)), 'a', newline='')
        else:
            f_ecif = open(os.path.join(tmpdata_dir, 'ecif_data', 'ECIF_{}.csv'.format(d)), 'a', newline='')
            f_elements = open(os.path.join(tmpdata_dir, 'ecif_data', 'ELEMENTS_{}.csv'.format(d)), 'a', newline='')
        # write csv headers
        possible_ecif = helper.get_possible_ecif()
        possible_elements = helper.get_possible_elements()
        ecif_header = ['PDB'] + possible_ecif
        elements_header = ['PDB'] + possible_elements
        ecif_writer = csv.writer(f_ecif)
        elements_writer = csv.writer(f_elements)
        ecif_writer.writerow(ecif_header)
        elements_writer.writerow(elements_header)
        for id in tqdm(ecif_ids):
            if id in ecif_core_ids:
                protein_file = os.path.join(ecif_core, id, "{}_protein.pdb".format(id))
                ligand_file = os.path.join(ecif_core, id, "{}_ligand.sdf".format(id))
                ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
                elements_list = helper.get_elements(protein_file, ligand_file, float(d))
                ecif_writer.writerow([id] + ecif_list)
                elements_writer.writerow([id] + elements_list)
                continue
            if id in ecif_2016_refined_ids:
                protein_file = os.path.join(ecif_2016_refined, id, "{}_protein.pdb".format(id))
                ligand_file = os.path.join(ecif_2016_refined, id, "{}_ligand.sdf".format(id))
                ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
                elements_list = helper.get_elements(protein_file, ligand_file, float(d))
                ecif_writer.writerow([id] + ecif_list)
                elements_writer.writerow([id] + elements_list)
                continue
            if id in ecif_2019_refined_ids:
                protein_file = os.path.join(ecif_2019_refined, id, "{}_protein.pdb".format(id))
                ligand_file = os.path.join(ecif_2019_refined, id, "{}_ligand.sdf".format(id))
                ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
                elements_list = helper.get_elements(protein_file, ligand_file, float(d))
                ecif_writer.writerow([id] + ecif_list)
                elements_writer.writerow([id] + elements_list)
                continue
            if id in ecif_2019_general_minus_refined_ids:
                protein_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_protein.pdb".format(id))
                ligand_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_ligand.sdf".format(id))
                ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
                elements_list = helper.get_elements(protein_file, ligand_file, float(d))
                ecif_writer.writerow([id] + ecif_list)
                elements_writer.writerow([id] + elements_list)
                continue
            not_listed_ids.add(id)
        f_ecif.close()
        f_elements.close()

    print("\n- Finished -\n")
    if len(not_listed_ids) > 0:
        print("Following IDs not listed in ECIF:\n")
        print(not_listed_ids)
        print("\n")


def collect_ecif_ex(distance_cutoffs, indices=None, tag=None):
    pass


def collect_ld(indices=None, tag=None):
    """

    :param indices:
    :param tag:
    :return:
    """
    if indices is not None:
        ecif_ids = indices
    else:
        ecif_dict = get_ecif_dict()
        ecif_ids = ecif_dict.keys()
    ecif_core_ids = IndexMng.get_index_from_dir(ecif_core)
    ecif_2016_refined_ids = IndexMng.get_index_from_dir(ecif_2016_refined)
    ecif_2019_refined_ids = IndexMng.get_index_from_dir(ecif_2019_refined)
    ecif_2019_general_minus_refined_ids = IndexMng.get_index_from_dir(ecif_2019_general_minus_refined)
    helper = ECIF(2016)
    lid_fts = {}
    not_listed_ids = set()
    print("\nCollecting Ligand Descriptors: \n")
    for id in tqdm(ecif_ids):
        if id in ecif_core_ids:
            ligand_file = os.path.join(ecif_core, id, "{}_ligand.sdf".format(id))
        elif id in ecif_2016_refined_ids:
            ligand_file = os.path.join(ecif_2016_refined, id, "{}_ligand.sdf".format(id))
        elif id in ecif_2019_general_minus_refined_ids:
            ligand_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_ligand.sdf".format(id))
        elif id in ecif_2019_refined_ids:
            ligand_file = os.path.join(ecif_2019_refined, id, "{}_ligand.sdf".format(id))
        else:
            not_listed_ids.add(id)
            continue
        if id not in lid_fts.keys():
            lid_fts[id] = helper.get_ligand_features_by_file(ligand_file)
    if tag is not None:
        f_ld = open(os.path.join(tmpdata_dir, 'ecif_data', 'RDKit_Descriptors_{}.csv'.format(tag)), 'a', newline='')
    else:
        f_ld = open(os.path.join(tmpdata_dir, 'ecif_data', 'RDKit_Descriptors.csv'), 'a', newline='')
    # write csv header
    ld_header = ['PDB'] + helper.get_ligand_descriptors()
    ld_writer = csv.writer(f_ld)
    ld_writer.writerow(ld_header)
    for id in ecif_ids:
        if id not in not_listed_ids:
            lid_ft = [id] + list(lid_fts[id])
            ld_writer.writerow(lid_ft)
    f_ld.close()
    print("\n- Finished -\n")
    if len(not_listed_ids) > 0:
        print("Following IDs not listed in ECIF:\n")
        print(not_listed_ids)
        print("\n")


def query_pk_by_id(qid):
    return float(affinity_data[affinity_data['pdbid'] == qid]['-logKd/Ki'])


def write_binding_data(tag, train_num):
    """
    for extra experiments, new BindingData.csv that associates with new fingerprint file and ligand descriptor file is needed
    :param tag: tag is needed
    :param train_num: size of train set, if not given, all data are labeled as 'Train'
    :return:
    """
    if tag is not None:
        ecif_data = pd.read_csv(os.path.join(tmpdata_dir, 'ecif_data', 'ECIF_6.0_{}.csv'.format(tag)), comment='#')
        f_bd = open("../Preprocess/BindingData_{}.csv".format(tag), 'a', newline='')
    else:
        ecif_data = pd.read_csv(os.path.join(tmpdata_dir, 'ecif_data', 'ECIF_6.0.csv'.format(tag)), comment='#')
        f_bd = open("../Preprocess/BindingData.csv".format(tag), 'a', newline='')
    bd_writer = csv.writer(f_bd)
    headers = ['PDB', 'SET', 'pK']
    bd_writer.writerow(headers)
    count = train_num if train_num is not None else len(ecif_data['PDB'])
    for id in ecif_data['PDB']:
        count -= 1
        # omotcha: this extra experiment assumes no split on dataset, and use auto test/train split when training
        if count >= 0:
            bd_writer.writerow([id, 'Train', query_pk_by_id(id)])
        else:
            bd_writer.writerow([id, 'Test', query_pk_by_id(id)])
    f_bd.close()


def helper_get_test():
    distance_cutoffs = '6.0'
    helper = ECIF(2016)
    protein_file = os.path.join(tmpdata_dir, "1a0q_protein.pdb")
    ligand_file = os.path.join(tmpdata_dir, "1a0q_ligandCD1.sdf")
    ecif_list = helper.get_ecif(protein_file, ligand_file, float(distance_cutoffs))
    elements_list = helper.get_elements(protein_file, ligand_file, float(distance_cutoffs))
    print("\n")
    print(ecif_list)


def best_cutoff_collecting():
    collect_ecif(['6.0'])
    collect_ld()


def all_cutoff_collecting():
    collect_ecif(['4.0', '4.5', '5.0', '5.5', '6.0', '6.5', '7.0', '7.5', '8.0', '8.5',
                  '9.0', '9.5', '10.0', '10.5', '11.0', '11.5', '12.0', '12.5', '13.0', '13.5',
                  '14.0', '14.5', '15.0'])
    collect_ld()


def best_cutoff_collecting_ign(tag):
    ids = get_all_ign_index()
    collect_ecif(['6.0'], indices=ids, tag=tag)
    collect_ld(indices=ids, tag=tag)


def best_cutoff_collecting_ign_2013(tag):
    n_test = 95
    ids = get_all_ign_index_2013(n_test)
    collect_ecif(['6.0'], indices=ids, tag=tag)
    collect_ld(indices=ids, tag=tag)


def my_test():
    # step 1:
    # best_cutoff_collecting_ign('ign')
    best_cutoff_collecting_ign_2013('test')

    # step 2:
    # write_binding_data('ign', train_num=8040)
    # write_binding_data('ign2013', train_num=8215)
    return


if __name__ == '__main__':
    my_test()
