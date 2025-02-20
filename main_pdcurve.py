from pypetting_extra import Experiment
from configurations import AssayConfiguration
from setup_pdcurve import pdSetup
from workflow_pdcurve import pdWorkflow
from general_classes import PathManager
import pandas as pd
import numpy as np

## use pyenv 3.12.0

folder = "twofold1to1"
exp_name = "test2"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

drug_pairs = [["amoxicillin", "colistin"], ["amoxicillin", "chloramphenicol"]]

pm = PathManager(basepath="current")
# drugs = pd.read_excel(pm.file_path("drugs.xlsx", folder="notes"))

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
)

config = AssayConfiguration()
setup = pdSetup(config, experiment, drug_pairs)

workflow = pdWorkflow(setup)
workflow.grow_overnight()
workflow.prepare_exp_cultures()
workflow.fill_replicates()
workflow.treat_cultures()
workflow.luminescence_read_loop()
