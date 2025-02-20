from pypetting_extra import Experiment
from configurations import AssayConfiguration
from setup_mk_pi import PiSetup, liha, roma, lid2_pos, plate2_pos
from workflow_mk_pi import PiWorkflow
from general_classes import PathManager
import pandas as pd
import numpy as np

## use pyenv 3.12.0

folder = "twofold1to1"
exp_name = "test_18022025"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
)

config = AssayConfiguration()
setup = PiSetup(config, experiment)

protocol = experiment.setup_protocol()
workflow = PiWorkflow(setup)
for drug_reservoir in setup.drug_reservoirs:
    print(drug_reservoir.name)

    wl = experiment.setup_worklist(f"{drug_reservoir.name}.gwl", protocol=protocol)
    wl.add(liha.sterile_wash())
    wl.add(
        roma.move_plate(
            drug_reservoir,
            lid2_pos,
            new_lid_gridsite=plate2_pos,
            end_with_covered_plate=False,
        )
    )
    wl.add(workflow.fill_LB(drug_reservoir))
    wl.add(workflow.dilution_row_drug_reservoir(drug_reservoir))
    wl.add(roma.incubate(drug_reservoir))
