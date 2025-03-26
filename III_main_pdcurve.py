from pypetting_extra import Experiment
from configurations import AssayConfiguration
from III_setup_pdcurve import pdSetup
from III_workflow_pdcurve import pdWorkflow
from general_classes import PathManager


## use pyenv 3.12.0

folder = "twofold1to1"
exp_name = "run_1"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

combination_idx = [0]

pm = PathManager(basepath="current")

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
    default_folders=["notes_III", "worklists_III"],
    default_keys=["notes_III", "wl_III"],
)

config = AssayConfiguration()
setup = pdSetup(config, experiment, combination_idx)


workflow = pdWorkflow(setup)
workflow.grow_overnight()
workflow.prefill_assay_plates()
workflow.prepare_exp_cultures()


workflow.treat_cultures()


workflow.luminescence_read_loop()
