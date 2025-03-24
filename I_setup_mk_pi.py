from configurations import AssayConfiguration
from pypetting_extra import Experiment
import pandas as pd
from pypetting.labware import labwares
from pypetting import Labware
from pypetting_extra import default_worktable as worktable
import numpy as np


storex = worktable.incubator
shelf = worktable.carrier["Shelf 8x4Pos"]

tip_arr1 = np.array([False, True, True, True] + 4 * [False])
tip_arr2 = np.array(4 * [False] + [True, True, True, False])
tip_arr = tip_arr1 + tip_arr2

liha = worktable.liha

mca = worktable.mca
roma = worktable.roma
tilter = worktable.tilter
lid1_pos = worktable.carrier["MP 2Pos Fixed"].gridsite(0)
lid2_pos = worktable.carrier["MP 2Pos Fixed"].gridsite(1)
lid3_pos = shelf.gridsite(30)
rotated_site = tilter.gridsite(0)
plate1_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(0)
plate2_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(1)
plate3_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(2)
corning_8row = Labware("8 Row DeepWell Corning", 8, 1)
deepwell96 = Labware("96 DeepWell Greiner", 8, 12)
greiner96 = labwares["greiner96"]
dw12col = Labware("12Col Through", 8, 12)


# Pi = Reservoir plate number i
class PiSetup:
    def __init__(self, config: AssayConfiguration, experiment: Experiment):
        self.config = config
        experiment.clone_folder("evoscripts")
        experiment.clone_folder("cmd_scripts")
        experiment.clone_folder("reader_settings")
        experiment.clone_folder("notes_I")
        experiment.clone_folder("plate_files")
        experiment.clone_folder("eval_code")
        self.lumread = experiment.setup_measurement("lum_count_exp.xml", "lum_files")
        experiment.setup_pickolo_folder("img", "img_files", liha)
        self.location = experiment.setup_location_file(folderkey="notes_I")
        self.experiment = experiment

        self.drug_reservoirs = []
        for i, drug in enumerate(config.drugs):
            plate = self.define_drug_reservoir(i, drug)
            self.drug_reservoirs.append(plate)

        self.medium_trough = self.define_medium_trough()

    def define_drug_reservoir(
        self,
        position,
        name,
        labware=dw12col,
    ):
        plate_name = f"{name}_reservoir"
        return storex.cartridges[4].define_plate(plate_name, labware, position)

    def define_medium_trough(
        self,
        site=0,
        name="medium trough",
        carrier_name="MP 3Pos Deck",
        labware=dw12col,
    ):
        plate_name = f"{name}"
        plate = worktable.carrier[carrier_name].define_plate(plate_name, labware, site)
        self.location.add(plate)
        return plate

    def storex_locations(self):
        return storex.locations
