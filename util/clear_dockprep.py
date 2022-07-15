"""
platform: win
env: any
name: clear_dockprep.py
"""
import os
from config import *


def clear_pocket_mol2(fp):
    """

    :param fp:
    :return:
    """
    for name in os.listdir(fp):
        if len(name) == 4 and os.path.isfile(fp + "%s\\%s_pocket.mol2" % (name, name)):
            os.remove(fp + "%s\\%s_pocket.mol2" % (name, name))


def remove_not_preprocessed_dir(fp):
    """

    :param fp:
    :return:
    """
    count = 0
    for name in os.listdir(fp):
        if len(name) == 4 and not os.path.isfile(fp + "%s\\%s_pocket.mol2" % (name, name)):
            for file in os.listdir(os.path.join(fp, name)):
                os.remove(os.path.join(fp, name, file))
            os.removedirs(os.path.join(fp, name))
            count += 1
    print("{} dirs removed".format(count))


if __name__ == '__main__':
    data_dir = general_except_refined_dir
    # remove_not_preprocessed_dir(data_dir)
