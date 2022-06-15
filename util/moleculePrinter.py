from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from config import *
import os


def print_mol(mol, name="mol"):
    opts = DrawingOptions()
    opts.includeAtomNumbers = True
    opts.bondLineWidth = 2
    img = Draw.MolToImage(mol, options=opts)
    img_dir = os.path.join(project_dir, "img", name+".PNG")
    img.save(img_dir)


def test():
    pass


if __name__ == '__main__':
    test()
