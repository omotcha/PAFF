"""
platform: any
env: any
name: IndexMng.py
Manage dataset index
"""
from config import *
import os


def get_index_from_dir(d):
    """
    get index from directory
    :param d: directory
    :return: index list
    """
    ids = []
    files = os.listdir(d)
    for f in files:
        if len(f) == 4:
            ids.append(f)
    return ids


def get_index_from_file(index_dot_txt):
    """
    get index from index file
    :param index_dot_txt: index file
    :return: index list
    """
    ids = []
    with open(index_dot_txt, 'r') as f:
        for line in f.readlines():
            ids.append(line[0:4])
    return ids


def test():
    """
    test
    :return:
    """
    # print(len(get_index_from_dir(ecif_2019_general_minus_refined)))
    ign_test = os.path.join(ign_index_dir, "test.index.txt")
    print(get_index_from_file(ign_test))


if __name__ == '__main__':
    test()
