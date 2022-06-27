from config import *
from util import moleculePrinter
from oddt.utils import *


def get_index():
    """
    link the pdbbind dataset, get the index from refined dataset and core dataset
    :return:
    refined dataset index list, core dataset index list
    """
    refined_ids = list(refined_dataset.sets['refined'].keys())
    core_ids = list(core_dataset.sets['core'].keys())
    # refined_ids = [i for i in refined_ids if i not in core_ids]
    return refined_ids, core_ids


def test():
    pass


if __name__ == '__main__':
    test()
