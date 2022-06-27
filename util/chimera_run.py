"""
run.py
"""
import os
from config import *


def main(fp):
    pdb_list = [k for k in os.listdir(fp)
                if len(k) == 4 and not os.path.isfile(fp + "%s\\%s_pocket.mol2" % (k, k))]
    print(len(pdb_list))
    count = 0
    for name in pdb_list:
        if len(name) != 4:
            continue
        PDB_file = file_path + name + '\\' + name
        process_dot_py = os.path.join(tmpdata_dir, 'chimera_process.py')
        os.system("chimera --nogui {} {}".format(process_dot_py, name))
        count += 1
        print(count)


file_path = refined_dir
main(file_path)
