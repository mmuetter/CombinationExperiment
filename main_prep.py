from pypetting_extra import Experiment
from configurations import AssayConfiguration
from setup_drug_plates import DrugPlateSetup
from workflow_prep import PrepWorkflow
from general_classes import PathManager
import pandas as pd
import numpy as np
from functions import layout_antibiotic_plate

## use pyenv 3.12.0

folder = "twofold1to1"
exp_name = "prep_plates_04022025"
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
workflow.prefill_drug_reservoir()
workflow.dilution_row_drug_reservoir()
prefix = "combine_"

base_df = layout_antibiotic_plate(drugs, 0.5, 0.5)
mic_dict = dict(zip(drugs.drug, drugs.micRef))
[conc for conc in config.concentration_gradient]
cols = config.concentration_gradient
cols.sort(reverse=True)
names = [f"antibiotics{"".join(str(s).split(".")) }" for s in cols]
concentration_col_dict = dict(
    zip(names, np.array(range(len(config.concentration_gradient))) + 1)
)
for plate, name, rel_conc in zip(
    setup.antibiotic_plates, config.plate_names, config.concentration_gradient
):
    src_col = concentration_col_dict[name]
    workflow.combine_drugs(plate, 0.5, 0.5, prefix + name, src_col)
    df = base_df.copy()
    df["rel_conc"] = rel_conc
    experiment.save_csv(df, "setup_" + name + ".csv")

drugs["cmax_mic"] = 32
drugs["cmax"] = drugs.cmax_mic * drugs.micRef
drugs["V_stock_ul"] = drugs.cmax / (drugs.stock)
experiment.save_csv(drugs, "fill_table.csv")
