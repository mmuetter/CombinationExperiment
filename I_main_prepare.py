from pypetting_extra import Experiment
from configurations import AssayConfiguration
from worktable import liha, roma, lid2_pos, plate2_pos, storex
from I_setup_prepare import PiSetup
from I_workflow_prepare import PiWorkflow
from files import SetupFiles
import pandas as pd

## use pyenv 3.12.0
folder = "twofold1to1"
exp_name = "run_5"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
    default_folders=["notes_I", "worklists_I", "plate_files"],
    default_keys=["notes_I", "wl_I", "plate_files"],
)

config = AssayConfiguration()

subreservoir_plan = pd.read_excel("subreservoir_plan.xlsx", index_col=0)
files = SetupFiles(config, experiment)
files.mk_fill_table()
files.mk_reservoir_files()
files.mk_combined_drug_files()
labels_df = files.mk_labels(subreservoir_plan)
files.save_labels_pdf(labels_df, "labels.pdf", cols=3)


for drug_abbrv in subreservoir_plan.index:
    print(drug_abbrv)
    num_A_plates = (
        subreservoir_plan.loc[drug_abbrv, "A_sub min"]
        + subreservoir_plan.loc[drug_abbrv, "A_sub extra"]
    )
    num_B_pairs = (
        subreservoir_plan.loc[drug_abbrv, "Bi min"]
        + subreservoir_plan.loc[drug_abbrv, "Bi extra"]
    )
    setup = PiSetup(config, experiment, num_A_plates, num_B_pairs, drug_abbrv)
    setup.storex_locations(save=True)
    workflow = PiWorkflow(setup, drug_abbrv)
    workflow.mk_drug_reservoir()
    workflow.mk_A_subreservoirs()
    workflow.mk_combination_plates_B()


# test
setup_test = PiSetup(config, experiment, 1, 1, "test")
workflow = PiWorkflow(setup_test, "test")
workflow.mk_A_subreservoirs()
workflow.mk_combination_plates_B()


#######################################################
# double check
asp_vol = {}
disp_vol = {}

for col in range(1, 13):
    disp_vol[col] = config.pi_LB_fill_volumes[col]
    feed_col = config.pi_feeding_cols[col]
    print(feed_col)
    if feed_col:
        transfer_vol = (
            config.pi_transfer_volumes[col] if config.pi_transfer_volumes[col] else 0
        )
        asp_vol[feed_col] = asp_vol.get(feed_col, 0) + transfer_vol
        disp_vol[col] = disp_vol.get(col, 0) + transfer_vol

total_vol = {}
for col in range(1, 13):
    total_vol[col] = disp_vol.get(col, 0) - asp_vol.get(col, 0)

print("asp_vol", asp_vol)
print("disp_vol", disp_vol)
print("total_vol", total_vol)
