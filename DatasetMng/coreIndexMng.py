import csv
import os
from config import *

"""

    Initially there is no index file in core set, so i have to create it.
    Suppose that the core set is a pure subset of refined set, 
    then the index file of core set can be generated through that of refined set.
    
"""


def get_core_ids():
    ids = []
    files = os.listdir(core_dir)
    for f in files:
        ids.append(f)
    return ids


def create_core_index():
    core_ids = get_core_ids()
    refined_index_dir = os.path.join(refined_dir, "index")
    refined_index_file = os.path.join(refined_index_dir, "INDEX_refined_data.2016")
    core_index_dir = os.path.join(core_dir, "index")
    core_index_file = os.path.join(core_index_dir, "INDEX_core_data.2016")
    core_index = []
    with open(refined_index_file, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0][0:4] in core_ids:
                core_index.append(row)
    with open(core_index_file, "w") as f:
        for row in core_index:
            f.writelines(str(row[0])+'\n')
    return


def test():
    pass


if __name__ == '__main__':
    test()
