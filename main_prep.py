from pypetting_extra import Experiment
from configurations import AssayConfiguration
from setupDrugPlates import DrugPlateSetup
from prepWorkflow import PrepWorkflow
from general_classes import PathManager
import pandas as pd
import numpy as np
from functions import layout_antibiotic_plate

## use pyenv 3.12.0

folder = "twofold1:1"
exp_name = "prep_plates_03022025"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

pm = PathManager(basepath="current")
drugs = pd.read_excel(pm.file_path("drugs.xlsx", folder="notes"))

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
)

config = AssayConfiguration()
setup = DrugPlateSetup(config, experiment)
workflow = PrepWorkflow(setup)
workflow._prefill_drug_reservoir()
prefix = "combine_1:1_"

base_df = layout_antibiotic_plate(drugs, 0.5, 0.5)
mic_dict = dict(zip(drugs.drug, drugs.micRef))

for plate, name, rel_conc in zip(
    setup.antibiotic_plates, config.plate_names, config.concentration_gradient
):
    workflow.combine_drugs(plate, 0.5, 0.5, prefix + name)
    df = base_df.copy()
    df["rel_conc"] = rel_conc
    experiment.save_csv(df, "setup_" + name + ".csv")
