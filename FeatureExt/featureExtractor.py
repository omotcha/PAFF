import inspect

from DatasetMng.datasetMngr import *
from tqdm import *
from config import *
from rdkit import Chem
import numpy as np
import os


def get_atom_type_prop(mol):
    # Since python 3.6, Dicts have become ordered
    props = {}
    if uses_mk:
        atoms_type = np.array([element.upper() for element in mol.atomtype])
    else:
        atoms = mol.GetAtoms()
        inspect.getmodule(atoms[0])
        for atom in atoms:
            atom_type = atom.GetSymbol()
            formal_charge = atom.GetFormalCharge()
            print("type: " + atom_type + " charge: " + str(formal_charge))
            # props["hydrophobic"] = atom_type == "C"


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
    ##########
    for i in tqdm(core_ids):
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


def test_get_atom():
    test_core_id = '3ao4'
    prot_f = os.path.join(core_dir, test_core_id, test_core_id + "_protein.pdb")
    ligd_f = os.path.join(core_dir, test_core_id, test_core_id + "_ligand.mol2")
    protein = Chem.MolFromPDBFile(prot_f)
    ligand = Chem.MolFromMol2File(ligd_f)
    get_atom_type_prop(protein)


def test_with_mk():
    """
    moleculekit is used for comparing results of each functions, and is excluded from this repo online
    :return:
    """
    pass


if __name__ == '__main__':
    test_get_atom()
