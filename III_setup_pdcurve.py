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
trough = Labware("Trough 300ml MCA96", 8, 12)


class pdSetup:
    def __init__(
        self, config: AssayConfiguration, experiment: Experiment, combination_indices
    ):
        self.config = config
        self.experiment = experiment

        # Soft Setup
        experiment.clone_folder("evoscripts")
        experiment.clone_folder("cmd_scripts")
        experiment.clone_folder("reader_settings")
        experiment.clone_folder("notes_III")
        experiment.clone_folder("eval_code")
        experiment.clone_folder("lum_files")
        self.lumread = experiment.setup_measurement("lum_count_exp.xml", "lum_files")
        self.location = experiment.setup_location_file(folderkey="notes_III")
        self.location.add(mca.tips[0])

        # Define Hardware
        self.combination_indices = combination_indices
        self.combinations = {}
        for key, combination in config.combinations.items():
            if key in combination_indices:
                print(f"Combination {key} is included.")
                self.combinations.update({key: combination})

        ## 384 well plate to distribute overnight cultures
        self.helper_plate = storex.define_plate_next_free_site(
            "helper_plate", greiner384
        )

        self.antibiotic_plates = self.define_antibiotic_plates(self.combinations)
        self.assay_plates = self.define_assay_plate(self.combinations)
        self.medium = self.define_medium_trough()
        self.overnight_12col = self.define_overnight_plate()
        experiment.save_csv(
            self.storex_locations(), "storex_locations.csv", folder_key="notes_III"
        )

    def define_overnight_plate(
        self,
        plate_name="overnight strains",
        labware=dw12col,
        position=0,
    ):
        plate = storex.cartridges[4].define_plate(plate_name, labware, position)
        return plate

    def define_medium_trough(
        self,
        plate_name="LB",
        position=2,
        carrier_name="MP 3Pos Deck",
        labware=trough,
    ):
        trough = worktable.carrier[carrier_name].define_plate(
            plate_name, labware, position
        )
        self.location.add(trough)
        return trough

    def define_assay_plate(
        self,
        combinations: dict,
    ):
        assayplates = []
        for _, combination in combinations.items():
            # I
            label = f"{combination['a']}_{combination['b']}_I"
            plate = storex.define_plate_next_free_site(label, greiner384)
            combination.update({"assay_I": plate})
            assayplates.append(plate)
            # II
            label = f"{combination['a']}_{combination['b']}_II"
            plate = storex.define_plate_next_free_site(label, greiner384)
            combination.update({"assay_II": plate})
            assayplates.append(plate)
        return assayplates

    def define_antibiotic_plates(
        self,
        combinations: dict,
    ):
        site = 31
        antibioticplates = []
        for _, combination in combinations.items():
            # I
            label = f"antibiotics_{combination['a']}_{combination['b']}_I"
            plate = shelf.define_plate(label, greiner96, site, is_store_pos=True)
            combination.update({"antibiotics_I": plate})
            self.location.add(plate)
            site -= 1
            antibioticplates.append(plate)
            # II
            label = f"antibiotics_{combination['a']}_{combination['b']}_II"
            plate = shelf.define_plate(label, greiner96, site, is_store_pos=True)
            combination.update({"antibiotics_II": plate})
            self.location.add(plate)
            site -= 1
            antibioticplates.append(plate)
        return antibioticplates

    def storex_locations(self):
        return storex.locations()
