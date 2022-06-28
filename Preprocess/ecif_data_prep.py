from config import *
import os
import csv
from DatasetMng import IndexMng
from util.ECIF import *
from tqdm import tqdm


def get_ecif_ids():
    ids = []
    with open('BindingData.csv', 'r') as f:
        f.readline()
        reader = csv.reader(f)
        for row in reader:
            ids.append(row[0])
    return ids


def collect_ecif(distance_cutoffs):
    ecif_ids = get_ecif_ids()
    ecif_core_ids = IndexMng.get_index_from_dir(ecif_core)
    ecif_2016_refined_ids = IndexMng.get_index_from_dir(ecif_2016_refined)
    ecif_2019_refined_ids = IndexMng.get_index_from_dir(ecif_2019_refined)
    ecif_2019_general_minus_refined_ids = IndexMng.get_index_from_dir(ecif_2019_general_minus_refined)
    helper = ECIF(2016)
    print("\nCollecting ECIF and ELEMENTS data: ")
    for d in distance_cutoffs:
        print("\n distance_cutoff: {}\n".format(d))
        f_ecif = open(os.path.join(tmpdata_dir, 'ecif_data', 'ECIF_{}.csv'.format(d)), 'a', newline='')
        f_elements = open(os.path.join(tmpdata_dir, 'ecif_data', 'ELEMENTS_{}.csv'.format(d)), 'a', newline='')
        # write csv headers
        possible_ecif = helper.get_possible_ecif()
        possible_elements = helper.get_possible_elements()
        ecif_header = ['PDB'] + possible_ecif
        elements_header = ['PDB'] + possible_elements
        ecif_writer = csv.writer(f_ecif)
        elements_writer = csv.writer(f_elements)
        ecif_writer.writerow(ecif_header)
        elements_writer.writerow(elements_header)
        for id in tqdm(ecif_ids):
            if id in ecif_core_ids:
                protein_file = os.path.join(ecif_core, id, "{}_protein.pdb".format(id))
                ligand_file = os.path.join(ecif_core, id, "{}_ligand.sdf".format(id))
                ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
                elements_list = helper.get_elements(protein_file, ligand_file, float(d))
                ecif_writer.writerow([id] + ecif_list)
                elements_writer.writerow([id] + elements_list)
                continue
            if id in ecif_2016_refined_ids:
                protein_file = os.path.join(ecif_2016_refined, id, "{}_protein.pdb".format(id))
                ligand_file = os.path.join(ecif_2016_refined, id, "{}_ligand.sdf".format(id))
                ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
                elements_list = helper.get_elements(protein_file, ligand_file, float(d))
                ecif_writer.writerow([id] + ecif_list)
                elements_writer.writerow([id] + elements_list)
                continue
            if id in ecif_2019_refined_ids:
                protein_file = os.path.join(ecif_2019_refined, id, "{}_protein.pdb".format(id))
                ligand_file = os.path.join(ecif_2019_refined, id, "{}_ligand.sdf".format(id))
                ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
                elements_list = helper.get_elements(protein_file, ligand_file, float(d))
                ecif_writer.writerow([id] + ecif_list)
                elements_writer.writerow([id] + elements_list)
                continue
            if id in ecif_2019_general_minus_refined_ids:
                protein_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_protein.pdb".format(id))
                ligand_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_ligand.sdf".format(id))
                ecif_list = helper.get_ecif(protein_file, ligand_file, float(d))
                elements_list = helper.get_elements(protein_file, ligand_file, float(d))
                ecif_writer.writerow([id] + ecif_list)
                elements_writer.writerow([id] + elements_list)
                continue
        f_ecif.close()
        f_elements.close()

    print("\n- Finished -\n")


def collect_ld():
    ecif_ids = get_ecif_ids()
    ecif_core_ids = IndexMng.get_index_from_dir(ecif_core)
    ecif_2016_refined_ids = IndexMng.get_index_from_dir(ecif_2016_refined)
    ecif_2019_refined_ids = IndexMng.get_index_from_dir(ecif_2019_refined)
    ecif_2019_general_minus_refined_ids = IndexMng.get_index_from_dir(ecif_2019_general_minus_refined)
    helper = ECIF(2016)
    lid_fts = {}
    print("\nCollecting Ligand Descriptors: \n")
    for id in tqdm(ecif_ids):
        if id in ecif_core_ids:
            ligand_file = os.path.join(ecif_core, id, "{}_ligand.sdf".format(id))
        elif id in ecif_2016_refined_ids:
            ligand_file = os.path.join(ecif_2016_refined, id, "{}_ligand.sdf".format(id))
        elif id in ecif_2019_general_minus_refined_ids:
            ligand_file = os.path.join(ecif_2019_general_minus_refined, id, "{}_ligand.sdf".format(id))
        elif id in ecif_2019_refined_ids:
            ligand_file = os.path.join(ecif_2019_refined, id, "{}_ligand.sdf".format(id))
        if id not in lid_fts.keys():
            lid_fts[id] = helper.get_ligand_features_by_file(ligand_file)
    f_ld = open(os.path.join(tmpdata_dir, 'ecif_data', 'RDKit_Descriptors.csv'), 'a', newline='')
    # write csv header
    ld_header = ['PDB'] + helper.get_ligand_descriptors()
    ld_writer = csv.writer(f_ld)
    ld_writer.writerow(ld_header)
    for id in ecif_ids:
        lid_ft = [id] + list(lid_fts[id])
        ld_writer.writerow(lid_ft)
    f_ld.close()
    print("\n- Finished -\n")


def my_test():
    distance_cutoffs = '6.0'
    helper = ECIF(2016)
    protein_file = os.path.join(tmpdata_dir, "1a0q_protein.pdb")
    ligand_file = os.path.join(tmpdata_dir, "1a0q_ligandCD1.sdf")
    ecif_list = helper.get_ecif(protein_file, ligand_file, float(distance_cutoffs))
    elements_list = helper.get_elements(protein_file, ligand_file, float(distance_cutoffs))
    print("\n")
    print(ecif_list)


def best_cutoff_collecting():
    collect_ecif(['6.0'])
    collect_ld()


def all_cutoff_collecting():
    collect_ecif(['4.0', '4.5', '5.0', '5.5', '6.0', '6.5', '7.0', '7.5', '8.0', '8.5',
                  '9.0', '9.5', '10.0', '10.5', '11.0', '11.5', '12.0', '12.5', '13.0', '13.5',
                  '14.0', '14.5', '15.0'])
    collect_ld()


def test():
    my_test()


if __name__ == '__main__':
    best_cutoff_collecting()
