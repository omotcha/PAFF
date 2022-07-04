from config import *
import os
import re


class AffinityDataMngr:
    _aff_data_loc = ""
    _aff_data = {}

    def _get_aff_data_from_index_file(self, f_index):
        """
        store affinity data from index file to a dictionary
        make sure that the index file named like INDEX_XXXXXX_data.XXXX, which has -pK column
        :param f_index: the index file
        :return:    the #items stored in the dictionary
        """
        with open(f_index, "r") as f:
            lines = f.readlines()
            for line in lines:
                seg = re.split(' +', line)
                if seg[0] != '#':
                    self._aff_data[seg[0]] = seg[3]
        return len(self._aff_data)

    # callables

    def get_pk_by_index(self, i):
        return self._aff_data[i]

    def get_aff_data(self):
        """
        get affinity data stored by list
        :return:
        """
        ret = []
        for k in self._aff_data.keys():
            ret.append([k, self._aff_data[k]])
        return ret

    def write_aff_data(self, file_name):
        """
        write affinity data to local for further steps
        :param file_name: it should be like affinity_data.csv
        :return:
        """
        with open(os.path.join(self._aff_data_loc, file_name), "w") as f:
            f.write('pdbid,-logKd/Ki\n')
            for k in self._aff_data.keys():
                f.write('{},{}\n'.format(k, self._aff_data[k]))

    def testADM(self):
        print(self._aff_data.values())

    def __init__(self, aff_loc, src):
        if aff_loc is None:
            self._aff_data_loc = tmpdata_dir
        else:
            self._aff_data_loc = aff_loc
        self._get_aff_data_from_index_file(src)


def test():
    aff_mngr = AffinityDataMngr(
        os.path.join(project_dir, "affdata"),
        os.path.join(extra_2020_index, "INDEX_general_PL_data.2020"))
    aff_mngr.write_aff_data("affinity_data_2020.csv")


if __name__ == '__main__':
    test()
