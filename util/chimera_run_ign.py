"""
platform: win
env: chimera_env
name: chimera_run_ign.py
"""
import os
from config import *


def main(finput, fpdb, fsdf, sc):
    names = [k for k in os.listdir(finput)
             if len(k) == 4 and not os.path.isfile(fpdb + "{}\\{}_protein.pdb".format(k, k))]
    entries = os.listdir(fp)
    files = [os.path.join(fp, entry) for entry in entries]
    count = 0
    for name in files:
        process_dot_py = "chimera_process_ign.py"
        os.system("chimera --nogui {} {}".format(process_dot_py, name))
        count += 1
        print(count)


input_path = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\coreset"
pdb_path = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\example\\mol2_files"
sdf_path = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\example\\sdf_files"
script = "D:\\AlexYu\\PAFF\\util\\chimera_process.py"
main(input_path, pdb_path, sdf_path, script)
