"""
platform: win+linux cluster
env: pafnucy_env
name: pafnucy_batch_predict.py
step 1: run script by run_directly() or save it first by save_script()
step 2: run (on cluster) predict.py@pafnucy with the output hdf file generated from step 1
step 3: rank predictions from the output csv file generated from step 2 and add smiles by add_smiles()
"""
import os
import sys
import pandas as pd
from rdkit import Chem


def get_cmd(script, pocket, lig_folder, output_folder):
    """

    :param script:
    :param pocket:
    :param lig_folder:
    :param output_folder:
    :return:
    """
    cmd = "python {} -l ".format(script)
    for f in os.listdir(lig_folder):
        cmd += os.path.join(lig_folder, f) + " "
    cmd += "--ligand_format sdf -p {} -o {}".format(pocket, os.path.join(output_folder, "complexes.hdf"))
    return cmd


def run_directly():
    print(len(sys.argv))
    if len(sys.argv) != 5:
        print("Argument Error.")
        print("Arguments of get_cmd: ")
        print("[prepare.py script location][pocket.mol2 file][ligand folder][output folder]")
        return
    cmd = get_cmd(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    os.system(cmd)


def save_script():
    prepare = "[prepare.py script location]"
    pocket = "[pocket.mol2 file]"
    ligands = "[ligand folder]"
    output = "[output folder]"
    cmd = get_cmd(prepare, pocket, ligands, output)
    with open(os.path.join(output, "script"), "w") as fout:
        fout.write(cmd)


def get_smiles_from_sdf(sdf):
    mols = [mol for mol in Chem.SDMolSupplier(sdf)]
    smi = Chem.MolToSmiles(mols[0])
    print(smi)


def add_smiles(f_pred, lig_folder):
    def get_smiles_from_file(f):
        mols = [mol for mol in Chem.SDMolSupplier(os.path.join(lig_folder, "{}.sdf".format(f)))]
        smi = Chem.MolToSmiles(mols[0])
        return smi
    pred = pd.read_csv(f_pred).sort_values(by="prediction", ascending=False)[0:1000]
    pred["smiles"] = pred.apply(lambda x: get_smiles_from_file(x["name"].astype(int)), axis=1)
    pred.to_csv("[prediction_with_smiles]")


if __name__ == '__main__':
    # predictions = "[all_predicions csv file]"
    # add_smiles(predictions, "[ligand folder]")
    pass
