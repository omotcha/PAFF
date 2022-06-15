from DatasetMng.datasetMngr import *
from tqdm import *
from config import *
from rdkit import Chem
from util import moleculePrinter

import os


def get_feature(protein_f, ligand_f):
    """
    get the complex feature from protein and ligand
    :param protein_f:
    :param ligand_f:
    :return:
    """
    protein = Chem.MolFromPDBFile(protein_f)
    ligand = Chem.MolFromMol2File(ligand_f)
    # testing molecule printer
    # moleculePrinter.print_mol(protein, "protein")
    # moleculePrinter.print_mol(ligand, "ligand")
    return 0


def get_core_features(core_ids):
    pids = []
    feats = []
    ##########
    # for test
    test_core_ids = ['3ao4']

    for i in tqdm(test_core_ids):
        prot_f = os.path.join(core_dir, i, i+"_protein.pdb")
        ligd_f = os.path.join(core_dir, i, i+"_ligand.mol2")
        if not os.path.isfile(prot_f) or not os.path.isfile(ligd_f):
            continue
        # get feature
        feat = get_feature(prot_f, ligd_f)
        pids.append(i)
        feats.append(feat)


def get_refined_feature(refined_ids):
    pass


def test():
    refined_ids, core_ids = get_index()
    get_core_features(core_ids)


def test_with_mk():
    pass


if __name__ == '__main__':
    test()
