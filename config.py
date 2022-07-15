"""
platform: win
env: any
name: config.py
"""
import os

project_dir = "D:\\AlexYu\\PAFF"
datasets_dir = "D:\\AlexYu\\datasets\\dataset"

"""
    dataset arrangement:    <datasets_dir> => {folder: pdbbind20xx}
                            <pdbbind20xx> => {folder:   "coreset", 
                                                        "refined-set", 
                                                        "general-minus-refined" | "general-set-except-refined"}
"""


tmpdata_dir = os.path.join(project_dir, "tmpdata")
util_dir = os.path.join(project_dir, "util")
preprocess_dir = os.path.join(project_dir, "Preprocess")
test_dir = os.path.join(project_dir, "tests")

# base dataset(or for pafnucy): PDBBind 2016
base_dir = os.path.join(datasets_dir, "pdbbind2016")
refined_dir = os.path.join(base_dir, "refined-set")
core_dir = os.path.join(base_dir, "coreset")
general_except_refined_dir = os.path.join(base_dir, "general-set-except-refined")

# for ECIF
ecif_2016_refined = refined_dir
ecif_2019_refined = os.path.join(datasets_dir, "pdbbind2019\\refined-set")
ecif_2019_general_minus_refined = os.path.join(datasets_dir, "pdbbind2019\\general-minus-refined")
ecif_core = core_dir
ecif_model_dir = os.path.join(project_dir, 'model', 'ecif')
ecif_data_dir = os.path.join(tmpdata_dir, 'ecif_data')

# for autogluon
autogluon_model_dir = os.path.join(project_dir, 'train', 'AutogluonModels')
ecif_example = os.path.join(autogluon_model_dir, "ag-20220629_084956\\models\\LightGBMXT\\model.pkl")

# for IGN
ign_index_dir = os.path.join(tmpdata_dir, 'ign_index')
core2013_dir = os.path.join(datasets_dir, "pdbbind2013\\coreset")
ign_example = os.path.join(base_dir, "example")

# for extra experiments using pdbbind2020
extra_2020_base = os.path.join(datasets_dir, "pdbbind2020")
extra_2020_index = os.path.join(extra_2020_base, "index")
extra_2020_refined = os.path.join(extra_2020_base, "refined-set")
extra_2020_general_minus_refined = os.path.join(extra_2020_base, "v2020-other-PL")
