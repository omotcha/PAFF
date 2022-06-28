"""
run.py
"""
import os
from config import *


def main(fp, sc):
    pdb_list = [k for k in os.listdir(fp)
                if len(k) == 4 and not os.path.isfile(fp + "%s\\%s_pocket.mol2" % (k, k))]
    print(len(pdb_list))
    count = 0
    for name in pdb_list:
        if len(name) != 4:
            continue
        PDB_file = file_path + name + '\\' + name
        process_dot_py = os.path.join(tmpdata_dir, sc)
        os.system("chimera --nogui {} {}".format(process_dot_py, name))
        count += 1
        print(count)


file_path = refined_dir
script = 'chimera_process.py'
main(file_path, script)
