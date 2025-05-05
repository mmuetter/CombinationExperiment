from I_setup_mk_pi import liha
from dataclasses import dataclass
import numpy as np


# Pi = Plate_i
class PiWorkflow:
    def __init__(self, setup: dataclass):
        self.experiment = experiment = setup.experiment
        self.setup = setup
        self.config = setup.config
        self.protocol = experiment.setup_protocol(folder_key="notes_I")
        self.tips = 8 * [True]
        self.tip_arr = [False] + 6 * [True] + [False]
        self.column_mask96_all = np.array(self.tips)
        self.column_mask96 = self.tip_arr
        self.column_mask384 = np.array(2 * [False] + 6 * [True, False] + 2 * [False])
        self.column_mask96_rotated = np.array([False] + 6 * [True] + 5 * [False])

    def fill_LB(self, drug_reservoir):
        wl = []
        for target_col, fill_volume in self.config.pi_LB_fill_volumes.items():
            wl += liha.vol_transfer(
                self.setup.medium_trough,
                1,
                drug_reservoir,
                target_col,
                fill_volume / sum(self.tips),
                self.column_mask96_all,
                tip_array=self.tips,
                liquid_class="Minimal FD ZMAX",
            )
        return wl

    def dilution_row_drug_reservoir(self, drug_reservoir):
        wl = []
        for target_col in self.config.pi_feeding_cols:
            wl += self.dilution_step(target_col, drug_reservoir)
            wl += liha.sterile_wash()
        return wl

    def dilution_step(self, target_col, drug_reservoir):
        src_col = self.config.pi_feeding_cols[target_col]
        total_vol = self.config.pi_transfer_volumes[target_col]
        wl = []
        print(src_col, target_col)
        if total_vol:
            transfer_volume = total_vol / sum(self.tips)
            wl += liha.vol_transfer(
                drug_reservoir,
                src_col,
                drug_reservoir,
                target_col,
                transfer_volume,
                self.column_mask96_all,
                tip_array=self.tips,
                liquid_class="Minimal FD ZMAX",
            )
            wl += liha.mix(
                drug_reservoir,
                target_col,
                250,
                2,
                self.column_mask96_all,
                tip_array=self.tips,
            )
        return wl

    def setup_worklist(self, name):
        return self.experiment.setup_worklist(
            name, protocol=self.protocol, key="wl_I", folder="worklists_I"
        )
