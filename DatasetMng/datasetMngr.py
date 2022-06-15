from config import *


def test():
    # core_ids = list(my_dataset.sets['core'].keys())
    # print(len(core_ids))
    refined_ids = list(refined_dataset.sets['refined'].keys())
    core_ids = list(core_dataset.sets['core'].keys())
    print(len(refined_ids))
    print(len(core_ids))


if __name__ == '__main__':
    test()
