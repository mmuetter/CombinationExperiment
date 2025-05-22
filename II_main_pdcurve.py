from pypetting_extra import Experiment
from configurations import AssayConfiguration
from II_setup_pdcurve import pdSetup
from II_workflow_pdcurve import pdWorkflow
from general_classes import PathManager


## use pyenv 3.12.0
folder = "twofold1to1"
exp_name = "run_5"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder
combination_idx = 14

pm = PathManager(basepath="current")

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
    default_folders=["notes_II", "worklists_II"],
    default_keys=["notes_II", "wl_II"],
)

config = AssayConfiguration()
setup = pdSetup(config, experiment, combination_idx)

workflow = pdWorkflow(setup)


workflow.grow_overnight()
workflow.prefill_assay_plates()
workflow.prepare_exp_cultures()
workflow.combine_plates()
workflow.treat_cultures()
workflow.luminescence_read_loop()
