from worktable import (
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
from dataclasses import dataclass
import numpy as np


# Pi = Plate_i
class PiWorkflow:
    def __init__(self, setup: dataclass, drug_abbrv: str):
        self.experiment = experiment = setup.experiment
        self.setup = setup
        self.config = setup.config
        self.drug_abbrv = drug_abbrv
        self.protocol = experiment.setup_protocol(folder_key="notes_I")
        self.tips = 8 * [True]
        self.tip_arr = [False] + 6 * [True] + [False]
        self.column_mask96_all = np.array(self.tips)
        self.column_mask96 = self.tip_arr
        self.column_mask384 = np.array(2 * [False] + 6 * [True, False] + 2 * [False])
        self.column_mask_rotated_I = np.array(6 * [True] + 6 * [False])
        self.column_mask_rotated_II = np.array(6 * [False] + 6 * [True])

    def mk_drug_reservoir(self):
        wl = self.setup_worklist(f"{self.drug_abbrv}_mk_drug_reservoir.gwl")
        wl.add(liha.sterile_wash())
        wl.add(
            roma.move_plate(
                self.setup.drug_reservoir,
                lid1_pos,
                new_lid_gridsite=plate1_pos,
                end_with_covered_plate=False,
            )
        )
        wl.add(
            self.fill_H2O(self.setup.drug_reservoir),
            msg=f"fill_water_to_drug_reservoir",
        )
        wl.add(
            self.dilution_row_drug_reservoir(self.setup.drug_reservoir),
            msg=f"perform_dilution_row_{self.setup.drug_reservoir.name}",
        )

    def fill_H2O(self, drug_reservoir):
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
                liquid_class="LB FD ZMAX",
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
                liquid_class="LB FD ZMAX",
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

    def mk_A_subreservoirs(self):
        wl = self.setup_worklist(f"{self.drug_abbrv}_mk_A_subreservoirs.gwl")
        for i, subreservoir_A in enumerate(self.setup.A_subreservoirs):
            self.add_drug_a(
                subreservoir_A, wl, i == len(self.setup.A_subreservoirs) - 1
            )

    def add_drug_a(self, subreservoir_A, wl, return_tips=False):
        reservoir_Pi = self.setup.drug_reservoir
        if not mca.tips_mounted:
            wl.add(mca.get_tips())
        transfer_vol = self.config.drug_sub_reservoir_vol
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
                subreservoir_A,
                plate2_pos,
                new_lid_gridsite=lid2_pos,
                end_with_covered_plate=False,
            )
        )
        wl.add(mca.aspirate(reservoir_Pi, 1, 1, transfer_vol))
        wl.add(
            mca.dispense(subreservoir_A, 1, 1, transfer_vol, retract=True),
            msg=f"adding_drug_A_to_{subreservoir_A.name}",
        )
        wl.add(roma.incubate(subreservoir_A))
        if return_tips:
            wl.add(mca.return_tips())

    def mk_combination_plates_B(self):
        wl = self.setup_worklist(f"{self.drug_abbrv}_mk_combination_plates_B_I.gwl")
        wl.add(liha.sterile_wash())
        for Pab_I in self.setup.combination_plate_B_I:
            self.add_drug_b(
                Pab_I,
                self.column_mask_rotated_I,
                wl,
            )
        wl.add(liha.sterile_wash())

        wl_II = self.setup_worklist(f"{self.drug_abbrv}_mk_combination_plates_B_II.gwl")
        for Pab_II in self.setup.combination_plate_B_II:
            self.add_drug_b(
                Pab_II,
                self.column_mask_rotated_II,
                wl_II,
            )
        wl_II.add(liha.sterile_wash())

    def add_drug_b(self, Pab_B, column_mask_rotated, wl):
        reservoir_Pi = self.setup.drug_reservoir
        transfer_vol = self.config.drug_plate_vol / 2
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
                Pab_B,
                plate2_pos,
                new_lid_gridsite=lid2_pos,
                end_with_covered_plate=False,
            )
        )
        wl.add(roma.move_roma_to_home())
        wl.add(
            liha.fill_plate(
                reservoir_Pi,
                Pab_B,
                transfer_vol,
                column_mask_rotated,
                self.column_mask96,
                tip_array=self.tip_arr,
                end_col=12,
                liquid_class="LB FD ZMAX",
                reverse_order=True,
                allow_multidispense=False,
            ),
            msg=f"adding_drug_B_to_{Pab_B.name}",
        )
        wl.add(roma.incubate(Pab_B))

    def setup_worklist(self, name):
        return self.experiment.setup_worklist(
            name, protocol=self.protocol, key="wl_I", folder_name="worklists_I"
        )
