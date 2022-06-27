import os
from config import *


def check_processed(fp):
    pdb_list = [k for k in os.listdir(fp)
                if len(k) == 4 and os.path.isfile(fp + "%s\\%s_pocket.mol2" % (k, k))]
    return len(pdb_list)


def check_not_preprocessed(fp):
    pdb_list = [k for k in os.listdir(fp)
                if len(k) == 4 and not os.path.isfile(fp + "%s\\%s_pocket.mol2" % (k, k))]
    return len(pdb_list)


def check_preprocessed(fp):
    print("checking {}".format(fp))
    print("{} processed".format(check_processed(fp)))
    print("{} not processed".format(check_not_preprocessed(fp)))


def check_inclusion(cd, rd):
    core_list = [k for k in os.listdir(cd)
                 if len(k) == 4 and os.path.isfile(cd + "%s\\%s_pocket.mol2" % (k, k))]
    refined_list = [k for k in os.listdir(rd)
                    if len(k) == 4 and os.path.isfile(rd + "%s\\%s_pocket.mol2" % (k, k))]
    for core in core_list:
        if core not in refined_list:
            # omotcha: in my test, only 3ge7 in core but not in refined(with xxxx_pocket.mol2), I removed it manually
            print(core)


if __name__ == '__main__':
    # check_preprocessed(general_minus_refined_dir)
    # check_inclusion(core_dir, refined_dir)
    print(check_not_preprocessed(general_except_core_dir))
