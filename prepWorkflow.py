from setup_pdcurve import (
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
    AssaySetup,
)
import numpy as np


class PrepWorkflow:
    def __init__(self, setup: AssaySetup):
        self.experiment = experiment = setup.experiment
        self.setup = setup
        self.config = setup.config
        self.protocol = experiment.setup_protocol()
        self.tip_arr = [False] + 6 * [True] + [False]
        self.column_mask96 = self.tip_arr
        self.column_mask384 = np.array(2 * [False] + 6 * [True, False] + 2 * [False])
        self.column_mask_rotated = np.array(6 * [False, True] + [False, False])

    def _prefill_drug_reservoir(self):
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
        wl.add(
            liha.dilution_row(
                self.setup.drug_reservoir,
                self.column_mask96,
                self.config.overnight_final_volume * 2,
                start_col=1,
                stop_at_col=len(self.config.concentration_gradient) - 1,
                dilution_factor=2,
                tip_array=self.tip_arr,
            )
        )
        wl.add(liha.sterile_wash())

    def combine_drugs(self, plate, pA, pB, label):
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
                self.column_mask_rotated,
                liquid_class="Minimal FD ZMAX",
                tip_array=self.tip_arr,
                start_col=2,
                end_col=7,
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
                self.column_mask_rotated,
                liquid_class="Minimal FD ZMAX",
                tip_array=self.tip_arr,
                start_col=2,
                end_col=7,
            )
        )
        wl.add(roma.store(plate))
        wl.add(liha.sterile_wash())

    def setup_worklist(self, name):
        return self.experiment.setup_worklist(name, protocol=self.protocol)
