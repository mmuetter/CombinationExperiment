from pypetting_extra import Experiment
from configurations import AssayConfiguration
from setup_drug_plates import CombineDrugsSetup
from workflow_prep import PrepWorkflow
from general_classes import PathManager
import pandas as pd
import numpy as np
from functions import layout_antibiotic_plate
import itertools


## use pyenv 3.12.0

folder = "twofold1to1"
exp_name = "test_new"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

pm = PathManager(basepath="current")
drugs = pd.read_excel(pm.file_path("drugs.xlsx", folder="notes"))

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
)

config = AssayConfiguration()
setup = CombineDrugsSetup(config, experiment)


# Example vector y
y = [1, 2, 3, 4]  # Replace with your actual vector

# Generate unique index pairs (i, j) with i < j to avoid duplicates
index_pairs = list(itertools.combinations(range(len(config.drugs)), 2))

# Construct dictionary with index as key and (i, j) pairs
combinations = {
    idx: {
        "i": i,
        "j": j,
        "a": config.drugs[i],
        "b": config.drugs[j],
        "Pa": setup.drug_reservoirs[i],
        "Pb": setup.drug_reservoirs[j],
    }
    for idx, (i, j) in enumerate(index_pairs)
}

setup.define_antibiotic_plates(combinations)
# Print result
print(combinations)


workflow = PrepWorkflow(setup)
for i, combination in combinations.items():
    workflow.add_drug_a(i, combination)
    workflow.add_drug_b(i, combination)
