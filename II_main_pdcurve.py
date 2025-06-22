from pypetting_extra import Experiment
from configurations import AssayConfiguration
from II_setup_pdcurve import pdSetup
from II_workflow_pdcurve import pdWorkflow
from general_classes import PathManager


## use pyenv 3.12.0
folder = "twofold1to1"
exp_name = "run_6"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

pm = PathManager(basepath="current")

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
    default_folders=["notes_II", "worklists_II"],
    default_keys=["notes_II", "wl_II"],
)


config = AssayConfiguration()
for combination_idx in range(0, 15):
    print(combination_idx)
    setup = pdSetup(config, experiment, combination_idx)
    workflow = pdWorkflow(setup)
    workflow.dilute_second_overnight()
    workflow.prefill_assay_plates()
    workflow.fill_helper_plate()
    workflow.infect_assays()
    workflow.combine_plates()
    workflow.treat_cultures()
    workflow.luminescence_read_loop()


self = workflow
wl = self.setup_worklist(f"test_mca.gwl")
wl.add(mca.get_pintool())
roma.move_plate(
    self.setup.helper_plate,
    plate2_pos,
    new_lid_gridsite=lid2_pos,
    end_with_covered_plate=False,
)
assay_plate = self.setup.assay_plates[0]
roma.move_plate(
    assay_plate,
    plate1_pos,
    new_lid_gridsite=lid1_pos,
    end_with_covered_plate=False,
)
wl.add(mca.replicate_with_pintool(self.setup.helper_plate, assay_plate, dip_offset=10))
wl.add(mca.move(self.setup.helper_plate))
wl.add(mca.drop_pintool())
