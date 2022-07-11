"""
platform: win
env: chimera_env
name: chimera_process_ign.py
"""

import chimera
from DockPrep import prep, AddH
from WriteMol2 import writeMol2
import sys
from config import *
import os

# name = sys.argv[-1]
mol2_file = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\example\\mol2_files\\1a0t_ligand.mol2"
sdf_path = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\example\\sdf_files"

mol = chimera.openModels.open(mol2_file)
"""
# optional:addHFunc=AddH.hbondAddHydrogens
# optional:addHFunc=AddH.simpleAddHydrogens
"""
# seems that hydrogen has already been added to ligand
prep(mol, addHFunc=AddH.hbondAddHydrogens, addCharges=False)
writeMol2(mol, "{}\\1a0t_ligand.sdf".format(sdf_path))
