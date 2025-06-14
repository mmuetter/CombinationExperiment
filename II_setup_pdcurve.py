from configurations import AssayConfiguration
from pypetting_extra import Experiment
import numpy as np
from worktable import mca, storex, shelf, carrier
from labware import (
    greiner384,
    greiner96,
    trough300 as trough,
    dw12col,
)


class pdSetup:
    def __init__(
        self, config: AssayConfiguration, experiment: Experiment, combination_idx
    ):
        self.config = config
        self.experiment = experiment

        # Soft Setup
        experiment.clone_folder("evoscripts")
        experiment.clone_folder("cmd_scripts")
        experiment.clone_folder("reader_settings")
        experiment.clone_folder("notes_II")
        experiment.clone_folder("eval_code")
        experiment.clone_folder("lum_files")
        self.combination_idx = combination_idx
        self.lumread = experiment.setup_measurement("lum_count_exp.xml", "lum_files")
        self.location = experiment.setup_location_file(folderkey="notes_II")
        self.location.add(mca.tips[0])

        # Define Hardware
        self.combination = config.combinations[combination_idx]

        ## 384 well plate to distribute overnight cultures
        self.helper_plate = storex.define_plate("helper_plate", greiner384, 0, 0)

        self.antibiotic_plates = self.define_antibiotic_plates()
        self.assay_plates = self.define_assay_plate()
        self.medium = self.define_medium_trough()
        self.overnight_12col = self.define_overnight_plate()
        self.next_overnight_12col = self.define_overnight_plate(
            "next_overnight_strains", dw12col, position=1
        )
        experiment.save_csv(
            self.storex_locations(), "storex_locations.csv", folder_key="notes_II"
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
        LB_trough = carrier[carrier_name].define_plate(plate_name, labware, position)
        self.location.add(LB_trough)
        return LB_trough

    def define_assay_plate(
        self,
    ):
        combination = self.combination
        assayplates = []

        # I
        label = f"{combination['a']}_{combination['b']}_I"
        plate = storex.define_plate(label, greiner384, 0, 1)
        combination.update({"assay_I": plate})
        assayplates.append(plate)
        # II
        label = f"{combination['a']}_{combination['b']}_II"
        plate = storex.define_plate(label, greiner384, 0, 2)
        combination.update({"assay_II": plate})
        assayplates.append(plate)
        return assayplates

    def define_antibiotic_plates(
        self,
    ):
        site = 31
        antibioticplates = []
        combination = self.combination

        # A_subreservoir
        label = f"antibiotics_{combination['a']}_{combination['b']}_A"
        plate = shelf.define_plate(label, greiner96, site, is_store_pos=True)
        combination.update({"subreservoir_A": plate})
        self.location.add(plate)
        site -= 1
        self.subreservoir_A = plate
        # B_combination_I
        label = f"antibiotics_{combination['a']}_{combination['b']}_I"
        plate = shelf.define_plate(label, greiner96, site, is_store_pos=True)
        combination.update({"antibiotics_I": plate})
        self.location.add(plate)
        site -= 1
        antibioticplates.append(plate)
        # B_combination_II
        label = f"antibiotics_{combination['a']}_{combination['b']}_II"
        plate = shelf.define_plate(label, greiner96, site, is_store_pos=True)
        combination.update({"antibiotics_II": plate})
        self.location.add(plate)
        antibioticplates.append(plate)
        site -= 1

        return antibioticplates

    def storex_locations(self):
        return storex.locations()
