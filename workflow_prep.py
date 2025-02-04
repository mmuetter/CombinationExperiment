from setup_drug_plates import (
    liha,
    roma,
    tilter,
    plate1_pos,
    plate2_pos,
    plate3_pos,
    lid1_pos,
    lid2_pos,
    lid3_pos,
    rotated_site,
)
from dataclasses import dataclass
import numpy as np


class PrepWorkflow:
    def __init__(self, setup: dataclass):
        self.experiment = experiment = setup.experiment
        self.setup = setup
        self.config = setup.config
        self.protocol = experiment.setup_protocol()
        self.tip_arr = [False] + 6 * [True] + [False]
        self.column_mask96 = self.tip_arr
        self.column_mask384 = np.array(2 * [False] + 6 * [True, False] + 2 * [False])
        self.column_mask96_rotated = np.array([False] + 6 * [True] + 5 * [False])

    def prefill_drug_reservoir(self):
        wl = self.setup_worklist("prefill_drug_reservoir.gwl")
        wl.add(liha.sterile_wash())

        wl.add(
            liha.fill_plate(
                self.setup.medium_reservoir,
                self.setup.drug_reservoir,
                self.config.antibiotic_reservoir_plate_final_vol,
                self.column_mask96,
                self.column_mask96,
                liquid_class="Minimal FD ZMAX",
                tip_array=self.tip_arr,
                start_col=2,
                end_col=len(self.config.concentration_gradient),
            )
        )

    def dilution_row_drug_reservoir(self):
        wl = self.setup_worklist("dilution_row_drug_reservoir.gwl")
        wl.add(
            liha.dilution_row(
                self.setup.drug_reservoir,
                self.column_mask96,
                self.config.antibiotic_reservoir_plate_final_vol * 2,
                start_col=1,
                stop_at_col=len(self.config.concentration_gradient) - 1,
                dilution_factor=2,
                tip_array=self.tip_arr,
            )
        )
        wl.add(liha.sterile_wash())

    def combine_drugs(self, plate, pA, pB, label, src_col):
        if pA + pB != 1:
            raise ValueError("pA + pB must be 1")
        wl = self.setup_worklist(label + ".gwl")
        wl.add(liha.sterile_wash())
        wl.add(
            roma.move_plate(
                plate,
                plate1_pos,
                new_lid_gridsite=lid1_pos,
                end_with_covered_plate=False,
            )
        )

        vol_A = pA * self.config.drug_plate_vol
        wl.add(
            liha.fill_plate(
                self.setup.drug_reservoir,
                plate,
                vol_A,
                self.column_mask96,
                self.column_mask96,
                liquid_class="Minimal FD ZMAX",
                tip_array=self.tip_arr,
                start_col=2,
                end_col=7,
                src_col=src_col,
            )
        )

        wl.add(
            roma.move_plate(
                plate,
                rotated_site,
                new_lid_gridsite=lid1_pos,
                end_with_covered_plate=False,
            )
        )
        vol_B = pB * self.config.drug_plate_vol
        wl.add(
            liha.fill_plate(
                self.setup.drug_reservoir,
                plate,
                vol_B,
                self.column_mask96,
                self.column_mask96_rotated,
                liquid_class="Minimal FD ZMAX",
                tip_array=self.tip_arr,
                start_col=2,
                end_col=7,
                src_col=src_col,
            )
        )
        wl.add(roma.store(plate))
        wl.add(liha.sterile_wash())

    def setup_worklist(self, name):
        return self.experiment.setup_worklist(name, protocol=self.protocol)
