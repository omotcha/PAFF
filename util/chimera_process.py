"""
process.py
"""

import chimera
from DockPrep import prep, AddH
from WriteMol2 import writeMol2
import sys
from config import *
import os

name = sys.argv[-1]
print(name)
file_path = refined_dir
PDB_file = file_path + name + '\\' + name
protein = chimera.openModels.open('%s_pocket.pdb' % PDB_file)
"""
# optional:addHFunc=AddH.hbondAddHydrogens
# optional:addHFunc=AddH.simpleAddHydrogens
"""
prep(protein, addHFunc=AddH.hbondAddHydrogens, addCharges=True)
writeMol2(protein, "%s_pocket.mol2" % PDB_file)
