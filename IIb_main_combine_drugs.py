from pypetting_extra import Experiment
from configurations import AssayConfiguration
from IIb_setup_combine import CombineDrugsSetup, storex
from II_workflow_combine import PrepWorkflow
from general_classes import PathManager
from II_worktable import liha, roma

## use pyenv 3.12.0
folder = "twofold1to1"
exp_name = "run_3"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

pm = PathManager(basepath="current")
# drugs = pd.read_excel(pm.file_path("drugs.xlsx", folder="notes"))

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
    default_folders=["notes_II", "worklists_II"],
    default_keys=["notes_II", "wl_II"],
)

config = AssayConfiguration()
combination_idx = 6
setup = CombineDrugsSetup(config, experiment, combination_idx)
combination = config.combinations[combination_idx]
setup.define_antibiotic_plate(combination_idx)
experiment.save_csv(storex.locations(), "storex_locations.csv", folder_key="notes_II")

workflow = PrepWorkflow(setup)
wl = workflow.setup_worklist(f"mk_combiplate_{combination_idx}.gwl")

drug_a = combination["a"]
reservoir_a = setup.drug_reservoirs[0]
# workflow.add_drug_a(reservoir_a, combination, wl)
# wl.add(roma.incubate(reservoir_a))

drug_b = combination["b"]
reservoir_b = setup.drug_reservoirs[1]
workflow.add_drug_b(reservoir_b, combination, wl)
wl.add(roma.incubate(reservoir_b))
wl.add(liha.sterile_wash())
