from IIb_setup_combine import (
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
        self.protocol = experiment.setup_protocol(folder_key="notes_II")
        self.tip_arr = np.array([False] + 6 * [True] + [False])
        self.column_mask96 = self.tip_arr
        self.column_mask384 = np.array(2 * [False] + 6 * [True, False] + 2 * [False])
        self.column_mask_rotated_I = np.array(6 * [True] + 6 * [False])
        self.column_mask_rotated_II = np.array(6 * [False] + 6 * [True])

    def add_drug_i(self, drug_i, combination_idx="all"):
        drug = self.config.drugs[drug_i]
        reservoir_Pi = self.setup.drug_reservoirs[drug_i]
        wl = self.setup_worklist(f"add_drug_{drug_i}.gwl")

        if combination_idx != "all":
            combinations = [self.config.combinations[combination_idx]]
        else:
            combinations = self.config.combinations.values()

        for combination in combinations:
            if drug == combination["a"]:
                print(f"Adding {drug} as component A")
                self.add_drug_a(reservoir_Pi, combination, wl)
            elif drug == combination["b"]:
                print(f"Adding {drug} as component B")
                self.add_drug_b(reservoir_Pi, combination, wl)

        wl.add(liha.sterile_wash())
        wl.add(roma.incubate(reservoir_Pi))
        wl.add(user_prompt("replace tips"))

    def add_drug_a(self, reservoir_Pi, combination, wl):
        wl.add(mca.get_tips())
        transfer_vol = self.config.drug_plate_vol / 2
        for Pab in [combination["Pab_I"], combination["Pab_II"]]:
            if reservoir_Pi.gridsite != lid1_pos:
                wl.add(
                    roma.move_plate(
                        reservoir_Pi,
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

            wl.add(mca.aspirate(reservoir_Pi, 1, 1, transfer_vol))
            wl.add(
                mca.dispense(Pab, 1, 1, transfer_vol),
                msg=f"adding_drug_A_to_{Pab.name}",
            )
            wl.add(roma.incubate(Pab))
        wl.add(mca.return_tips())

    def add_drug_b(self, reservoir_Pi, combination, wl):
        transfer_vol = self.config.drug_plate_vol / 2
        wl.add(liha.sterile_wash())
        for Pab, src_col_mask in zip(
            [combination["Pab_II"], combination["Pab_I"]],
            [self.column_mask_rotated_II, self.column_mask_rotated_I],
        ):
            if reservoir_Pi.gridsite != rotated_site:
                wl.add(
                    roma.move_plate(
                        reservoir_Pi,
                        rotated_site,
                        new_lid_gridsite=plate1_pos,
                        end_with_covered_plate=False,
                        dest_rotated=True,
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

            wl.add(roma.move_roma_to_home())

            wl.add(
                liha.fill_plate(
                    reservoir_Pi,
                    Pab,
                    transfer_vol,
                    src_col_mask,
                    self.column_mask96,
                    tip_array=self.tip_arr,
                    end_col=12,
                    liquid_class="Minimal FD ZMAX",
                ),
                msg=f"adding_drug_B_to_{Pab.name}",
            )
            wl.add(roma.incubate(Pab))

    def setup_worklist(
        self,
        name,
    ):
        return self.experiment.setup_worklist(
            name, protocol=self.protocol, key="wl_II", folder_name="worklists_II"
        )
