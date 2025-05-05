from pypetting_extra import Experiment
from configurations import AssayConfiguration
from I_setup_mk_pi import PiSetup, liha, roma, lid2_pos, plate2_pos, storex
from I_workflow_mk_pi import PiWorkflow
from files import SetupFiles


## use pyenv 3.12.0

folder = "twofold1to1"
exp_name = "run_3"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
    default_folders=["notes_I", "worklists_I"],
    default_keys=["notes_I", "wl_I"],
)

config = AssayConfiguration()

setup = PiSetup(config, experiment)
files = SetupFiles(config, experiment)
files.mk_fill_table()
files.mk_reservoir_files()
files.mk_combined_drug_files()

experiment.save_csv(storex.locations(), "storex_locations.csv", folder_key="notes_I")

protocol = experiment.setup_protocol(folder_key="notes_I")
workflow = PiWorkflow(setup)
for i, drug_reservoir in enumerate(setup.drug_reservoirs):
    print(drug_reservoir.name)

    wl = experiment.setup_worklist(
        f"reservoir_{i}.gwl", protocol=protocol, key="wl_I", folder_name="worklists_I"
    )
    wl.add(liha.sterile_wash())
    wl.add(
        roma.move_plate(
            drug_reservoir,
            lid2_pos,
            new_lid_gridsite=plate2_pos,
            end_with_covered_plate=False,
        )
    )
    wl.add(workflow.fill_LB(drug_reservoir), msg=f"fill_LB_plate_{drug_reservoir.name}")
    wl.add(
        workflow.dilution_row_drug_reservoir(drug_reservoir),
        msg=f"perform_dilution_row_{drug_reservoir.name}",
    )
    wl.add(roma.incubate(drug_reservoir), msg="dilution_row_done")
