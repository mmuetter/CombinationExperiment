from configurations import AssayConfiguration
from pypetting_extra import Experiment
from worktable import (
    liha,
    storex,
    shelf,
    worktable,
)
from labware import rotated_12col, greiner96, dw12col


# Pi = Reservoir plate number i
class PiSetup:
    def __init__(
        self,
        config: AssayConfiguration,
        experiment: Experiment,
        num_A_plates: int,
        num_B_pairs: int,
        drug_abbrv: str,
    ):
        self.config = config
        experiment.clone_folder("evoscripts")
        experiment.clone_folder("cmd_scripts")
        experiment.clone_folder("reader_settings")
        experiment.clone_folder("notes_I")
        experiment.clone_folder("plate_files")
        experiment.clone_folder("eval_code")
        self.lumread = experiment.setup_measurement("lum_count_exp.xml", "lum_files")
        experiment.setup_pickolo_folder("img", "img_files", liha)
        self.location = experiment.setup_location_file(
            folderkey="notes_I", filename=f"{drug_abbrv}_locations.csv"
        )
        self.experiment = experiment

        self.drug_reservoir = self.define_drug_reservoir()
        self.medium_trough = self.define_medium_trough()
        self.site = 31
        self.A_subreservoirs = []
        self.combination_plate_B_I = []
        self.combination_plate_B_II = []
        for site_idx in range(num_A_plates):
            self.A_subreservoirs.append(self.define_subreservoir_A(site_idx))

        for site_idx in range(num_B_pairs):
            plate_I, plate_II = self.define_antibiotic_combination_plates_B(site_idx)
            self.combination_plate_B_I.append(plate_I)
            self.combination_plate_B_II.append(plate_II)

    def define_drug_reservoir(
        self,
        position=0,
        name="drug_reservoir",
        labware=dw12col,
    ):
        plate_name = f"{name}_reservoir"
        reservoir = storex.cartridges[4].define_plate(plate_name, labware, position)
        reservoir.add_rotated_labw(rotated_12col)
        return reservoir

    def define_antibiotic_combination_plates_B(
        self,
        site_idx,
        name="antibiotic_combination_plate",
        labware=greiner96,
        I_cart_idx=1,
        II_cart_idx=2,
    ):
        plate_I = storex.define_plate(f"{name}_I", labware, I_cart_idx, site_idx)
        plate_II = storex.define_plate(f"{name}_II", labware, II_cart_idx, site_idx)
        return plate_I, plate_II

    def define_subreservoir_A(
        self, site_idx, name="subreservoir_A", labware=greiner96, cart_idx=0
    ):
        sub_reservoir = storex.define_plate(name, labware, cart_idx, site_idx)
        return sub_reservoir

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

    def storex_locations(self, save=False):
        if save:
            self.experiment.save_csv(
                storex.locations(), "storex_locations.csv", folder_key="notes_I"
            )
        return storex.locations()
