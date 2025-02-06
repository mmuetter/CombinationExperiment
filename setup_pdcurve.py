from configurations import AssayConfiguration
from pypetting_extra import Experiment
import re
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
rotated_site = tilter.gridsite(0)
plate1_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(0)
plate2_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(1)
plate3_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(2)
corning_8row = Labware("8 Row DeepWell Corning", 8, 1)
greiner96 = labwares["greiner96"]
greiner384 = labwares["greiner384"]
dw12col = Labware("12Col Through", 8, 12)


class pdSetup:
    def __init__(
        self, config: AssayConfiguration, experiment: Experiment, antibiotic_plate_names
    ):
        self.config = config
        self.experiment = experiment

        # Soft Setup
        experiment.clone_folder("evoscripts")
        experiment.clone_folder("cmd_scripts")
        experiment.clone_folder("reader_settings")
        experiment.clone_folder("notes")
        experiment.clone_folder("eval_code")
        experiment.clone_folder("lum_files")
        self.lumread = experiment.setup_measurement("lum_count_exp.xml", "lum_files")
        self.location = experiment.setup_location_file()
        self.location.add(mca.tips[0])

        # Define Hardware
        self.antibiotic_plate_names = antibiotic_plate_names
        self.assay_plate_names = [
            f"assay{re.search(r'\d+', s).group()}" for s in antibiotic_plate_names
        ]
        self.overnight_12col = self.define_overnight_plate()
        self.medium_trough = self.define_medium_trough()
        self.assayplates = self.define_assay_plate(self.assay_plate_names)
        self.antibiotic_plates = self.define_antibiotic_plates(antibiotic_plate_names)
        self.assay_dict = dict(zip(self.assay_plate_names, self.assayplates))
        experiment.save_csv(self.storex_locations(), "storex_locations.csv")

    def define_overnight_plate(
        self,
        plate_name="overnight strains",
        labware=dw12col,
        position=0,
    ):
        return storex.cartridges[4].define_plate(plate_name, labware, position)

    def define_medium_trough(
        self,
        plate_name="LB",
        position=2,
        carrier_name="MP 3Pos Deck",
        labware=corning_8row,
    ):
        trough = worktable.carrier[carrier_name].define_plate(
            plate_name, labware, position
        )
        self.location.add(trough)
        return trough

    def define_antibiotic_plates(self, plate_names):
        plates = []
        for i, name in enumerate(plate_names):
            print(i, name)
            plate = shelf.define_plate(name, greiner96, 31 - i, is_store_pos=True)
            self.location.add(plate)
            plates.append(plate)
        return plates

    def define_assay_plate(self, plate_names):
        plates = []
        for i, name in enumerate(plate_names):
            plate = storex.define_plate(name, greiner384, 0, i)
            plate.store_pos = shelf.gridsite(31 - len(plate_names) - i)
            self.location.add(plate)
            plates.append(plate)
        return plates

    def storex_locations(self):
        return storex.locations()
