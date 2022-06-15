from DatasetMng.datasetMngr import *
from tqdm import *
from config import *
import os


def get_feature(protein, ligand):
    """
    get the complex feature from protein and ligand
    :param protein:
    :param ligand:
    :return:
    """
    pass


def get_core_features(core_ids):
    pids = []
    feats = []
    for i in tqdm(core_ids):
        prot_f = os.path.join(core_dir, i, i+"_protein.pdb")
        ligd_f = os.path.join(core_dir, i, i+"_ligand.mol2")
        if not os.path.isfile(prot_f) or not os.path.isfile(ligd_f):
            continue
        #############
        # get feature

        #############
        pids.append(i)
        # should be replaced with feature
        feats.append(i)


def get_refined_feature(refined_ids):
    pass


def test():
    refined_ids, core_ids = get_index()
    get_core_features(core_ids)


if __name__ == '__main__':
    test()
