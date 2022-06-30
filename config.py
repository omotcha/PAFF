import os

project_dir = "D:\\AlexYu\\PAFF"
tmpdata_dir = "D:\\AlexYu\\PAFF\\tmpdata"
util_dir = "D:\\AlexYu\\PAFF\\util"
base_dir = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016"
refined_dir = os.path.join(base_dir, "refined-set")
core_dir = os.path.join(base_dir, "coreset")
general_except_core_dir = os.path.join(base_dir, "general-set-except-refined")


# for ECIF
ecif_2016_refined = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\refined-set"
ecif_2019_refined = "D:\\AlexYu\\datasets\\dataset\\pdbbind2019\\refined-set"
ecif_2019_general_minus_refined = "D:\\AlexYu\\datasets\\dataset\\pdbbind2019\\general-minus-refined"
ecif_core = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\coreset"
ecif_model_dir = os.path.join(project_dir, 'model', 'ecif')
ecif_data_dir = os.path.join(tmpdata_dir, 'ecif_data')

# for autogluon
autogluon_model_dir = os.path.join(project_dir, 'train', 'AutogluonModels')
ecif_example = "D:\\AlexYu\\PAFF\\train\\AutogluonModels\\ag-20220629_084956\\models\\LightGBMXT\\model.pkl"

# for IGN
ign_index_dir = os.path.join(tmpdata_dir, 'ign_index')
core2013_dir = "D:\\AlexYu\\datasets\\dataset\\pdbbind2013\\coreset"
