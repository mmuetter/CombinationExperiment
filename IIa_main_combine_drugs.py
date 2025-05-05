from pypetting_extra import Experiment
from configurations import AssayConfiguration
from II_setup_combine import CombineDrugsSetup, storex
from II_workflow_combine import PrepWorkflow
from general_classes import PathManager

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
setup = CombineDrugsSetup(config, experiment)
combinations = config.combinations

setup.define_antibiotic_plates(combinations)
experiment.save_csv(storex.locations(), "storex_locations.csv", folder_key="notes_II")

workflow = PrepWorkflow(setup)
for i, drug in enumerate(config.drugs):
    print("\n" + drug)
    workflow.add_drug_i(i, combinations)
