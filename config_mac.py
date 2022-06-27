import os
from oddt import datasets

project_dir = "/Users/o2cy/PycharmProjects/PAFF"

base_dir = "/Users/o2cy/Experiments/kdeep/datasets"
refined_dir = os.path.join(base_dir, "refined-set")
refined_dataset = datasets.pdbbind(home=refined_dir, default_set="refined", version=2016)

core_dir = os.path.join(base_dir, "core-set")
core_dataset = datasets.pdbbind(home=core_dir, default_set="core", version=2016)

uses_mk = False
