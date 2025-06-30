from pypetting_extra import Experiment
from configurations import AssayConfiguration
from II_setup_pdcurve import pdSetup
from II_workflow_pdcurve import pdWorkflow
from general_classes import PathManager


## use pyenv 3.12.0
folder = "twofold1to1"
exp_name = "run_7"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

pm = PathManager(basepath="current")

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
    default_folders=["notes_II", "worklists_II"],
    default_keys=["notes_II", "wl_II"],
)


## ATTENTION: This scripts contains exp-inc duration of 2.5 hours. Double check if you really want that.
# user prompt:
input(
    "This script contains an incubation time of 2.5 hours. Do you want to continue? (y/n): "
)

config = AssayConfiguration()
for combination_idx in range(0, 15):
    print(combination_idx)
    setup = pdSetup(config, experiment, combination_idx)
    workflow = pdWorkflow(setup)
    workflow.prefill_assay_plates()
    workflow.fill_helper_plate()
    workflow.infect_assays()
    workflow.mk_next_day_helper()
    workflow.combine_plates()
    workflow.treat_cultures()
    workflow.luminescence_read_loop()
