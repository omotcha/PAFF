from rdkit import Chem
import os


def read_ligand_from_dir(base_dir, id):
    # first attempt: read ligand.mol2
    ligand_file = os.path.join(base_dir, "{}_ligand.mol2".format(id))
    ligand = Chem.MolFromMol2File(ligand_file)
    if ligand is None:
        # second attempt: read ligand.sdf to smiles
        suppl = Chem.SDMolSupplier(os.path.join(base_dir, "{}_ligand.sdf".format(id)))
        smiles = [Chem.MolToSmiles(mol) for mol in suppl if mol]
        if len(smiles) != 0:
            # convert smiles back to rdkit.molecule, only the first smiles been converted
            ligand = Chem.MolFromSmiles(smiles[0])
            return ligand
        else:
            return None


def test():
    pass


if __name__ == '__main__':
    test()
