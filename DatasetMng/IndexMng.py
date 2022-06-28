from config import *
import os


def get_index_from_dir(d):
    ids = []
    files = os.listdir(d)
    for f in files:
        if len(f) == 4:
            ids.append(f)
    return ids


def test():
    print(len(get_index_from_dir(ecif_2019_general_minus_refined)))


if __name__ == '__main__':
    test()
