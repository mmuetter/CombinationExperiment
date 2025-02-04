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


class DrugPlateSetup:
    def __init__(self, config: AssayConfiguration, experiment: Experiment):
        self.config = config
        experiment.clone_folder("evoscripts")
        experiment.clone_folder("cmd_scripts")
        experiment.clone_folder("reader_settings")
        experiment.clone_folder("notes")
        experiment.clone_folder("eval_code")
        self.lumread = experiment.setup_measurement("lum_count_exp.xml", "lum_files")
        experiment.setup_pickolo_folder("img", "img_files", liha)
        self.location = experiment.setup_location_file()
        self.experiment = experiment
        self.drug_reservoir = self.define_reservoir_plate()
        self.medium_reservoir = self.define_medium_trough()
        self.antibiotic_plates = self.define_antibiotic_plates()

    def define_reservoir_plate(
        self,
        plate_name="antibiotic reservoir",
        labware=deepwell96,
        carrier_name="MP 3Pos Deck",
        position=0,
    ):
        plate = worktable.carrier[carrier_name].define_plate(
            plate_name, labware, position
        )
        self.location.add(plate)
        return plate

    def define_medium_trough(
        self,
        plate_name="medium trough",
        position=1,
        carrier_name="MP 3Pos Deck",
        labware=corning_8row,
    ):
        plate = worktable.carrier[carrier_name].define_plate(
            plate_name, labware, position
        )
        self.location.add(plate)
        return plate

    def define_antibiotic_plates(self):
        plates = []
        for i, label in enumerate(self.config.plate_names):
            plate = shelf.define_plate(label, greiner96, 31 - i, is_store_pos=True)
            plates.append(plate)
            self.location.add(plate)
        return plates

    def storex_locations(self):
        return storex.locations
