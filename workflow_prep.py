from setup_drug_plates import (
    liha,
    roma,
    mca,
    tilter,
    plate1_pos,
    plate2_pos,
    plate3_pos,
    lid1_pos,
    lid2_pos,
    lid3_pos,
    rotated_site,
)
from pypetting import user_prompt
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
        self.column_mask_rotated_I = np.array(6 * [True] + 6 * [False])
        self.column_mask_rotated_II = np.array(6 * [True] + 6 * [False])

    def add_drug_a(self, i, combination):
        wl = self.setup_worklist(f"combination_{i}_add_a.gwl")
        transfer_vol = self.config.drug_plate_vol / 2
        Pa = combination["Pa"]
        wl.add(mca.get_tips())
        for Pab in [combination["Pab_I"], combination["Pab_II"]]:

            wl.add(
                roma.move_plate(
                    Pa,
                    lid1_pos,
                    new_lid_gridsite=plate1_pos,
                    end_with_covered_plate=False,
                )
            )

            wl.add(
                roma.move_plate(
                    Pab,
                    plate2_pos,
                    new_lid_gridsite=lid2_pos,
                    end_with_covered_plate=False,
                )
            )

            wl.add(mca.aspirate(Pa, 1, 1, transfer_vol))
            wl.add(mca.dispense(Pab, 1, 1, transfer_vol))
            wl.add(roma.incubate(Pab))
        wl.add(mca.return_tips())
        wl.add(roma.incubate(Pa))
        wl.add(user_prompt("replace tips"))

    def add_drug_b(self, i, combination):
        wl = self.setup_worklist(f"combination_{i}_add_b.gwl")
        transfer_vol = self.config.drug_plate_vol / 2
        Pb = combination["Pb"]
        wl.add(liha.sterile_wash())
        for Pab, src_col_mask in zip(
            [combination["Pab_II"], combination["Pab_I"]],
            [self.column_mask_rotated_II, self.column_mask_rotated_I],
        ):
            wl.add(
                roma.move_plate(
                    Pb,
                    rotated_site,
                    new_lid_gridsite=lid3_pos,
                    end_with_covered_plate=False,
                )
            )
            wl.add(
                roma.move_plate(
                    Pab,
                    plate2_pos,
                    new_lid_gridsite=lid2_pos,
                    end_with_covered_plate=False,
                )
            )
            wl.add(
                liha.fill_plate(
                    Pb,
                    Pab,
                    transfer_vol,
                    src_col_mask,
                    self.column_mask96,
                    tip_array=self.tip_arr,
                    end_col=8,
                )
            )
            wl.add(roma.incubate(Pb))
        wl.add(liha.sterile_wash())
        wl.add(roma.incubate(Pb))

    def setup_worklist(self, name):
        return self.experiment.setup_worklist(name, protocol=self.protocol)
