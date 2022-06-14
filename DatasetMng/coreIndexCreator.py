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


def get_refined_dict():
    core_ids = get_core_ids()
    refined_index_dir = os.path.join(refined_dir, "index")
    refined_index_file = os.path.join(refined_index_dir, "INDEX_refined_data.2016")
    rows = []
    with open(refined_index_file, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            print(row[0])
            if row[0] in core_ids:
                rows.append(row)

    print(len(rows))
    return


def test():
    get_refined_dict()


if __name__ == '__main__':
    test()
