import os
from oddt import datasets

project_dir = "D:\\AlexYu\\PAFF"
tmpdata_dir = "D:\\AlexYu\\PAFF\\tmpdata"
util_dir = "D:\\AlexYu\\PAFF\\util"
base_dir = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016"
refined_dir = os.path.join(base_dir, "refined-set")
refined_dataset = datasets.pdbbind(home=refined_dir, default_set="refined", version=2016)

core_dir = os.path.join(base_dir, "coreset")
# core_dataset = datasets.pdbbind(home=core_dir, default_set="core", version=2016)

general_except_core_dir = os.path.join(base_dir, "general-set-except-refined")
