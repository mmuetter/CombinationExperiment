from worktable import (
    liha,
    roma,
    mca,
    plate1_pos,
    plate2_pos,
    plate3_pos,
    lid1_pos,
    lid2_pos,
    shaker_pos,
)
from dataclasses import dataclass
import numpy as np
from pypetting import start_timer, wait_timer, user_prompt
from icecream import ic
from itertools import product


class pdWorkflow:
    def __init__(self, setup: dataclass):
        self.experiment = experiment = setup.experiment
        self.setup = setup
        self.combination = setup.combination
        self.config = setup.config
        self.protocol = experiment.setup_protocol(
            file_name=f"timelog_c{setup.combination_idx}.csv", folder_key="notes_II"
        )
        self.combination_idx = setup.combination_idx
        self.tip_arr = [False] + 6 * [True] + [False]
        self.tip_arr8 = 8 * [True]
        self.column_mask96_8 = 8 * [True]
        self.column_mask96 = self.tip_arr
        self.column_mask384_A = np.array(2 * [False] + 6 * [True, False] + 2 * [False])
        self.column_mask384_B = np.array(2 * [False] + 6 * [False, True] + 2 * [False])

    def combine_plates(self):
        wl = self.setup_worklist(f"c{self.combination_idx}_combine_plates.gwl")
        plate_A = self.setup.subreservoir_A

        wl.add(
            roma.move_plate(
                plate_A,
                lid2_pos,
                new_lid_gridsite=plate2_pos,
                end_with_covered_plate=False,
            )
        )

        for _, suffix in zip([1, 0], ["_II", "_I"]):
            wl.add(mca.get_tips())
            plate_B = self.combination["antibiotics" + suffix]
            wl.add(
                roma.move_plate(
                    plate_B,
                    plate3_pos,
                    new_lid_gridsite=lid1_pos,
                    end_with_covered_plate=False,
                )
            )
            wl.add(
                mca.aspirate(
                    plate_A,
                    1,
                    1,
                    self.config.drug_plate_vol / 2,
                    liquid_class="Water free dispense",
                )
            )
            wl.add(
                mca.dispense(
                    plate_B,
                    1,
                    1,
                    self.config.drug_plate_vol / 2,
                    liquid_class="Water free dispense",
                ),
                msg="antibiotic_plates_combined",
            )
            wl.add(roma.store(plate_B))
            wl.add(mca.return_tips())

        wl.add(roma.store(plate_A))

    def fill_helper_plate(self):
        wl = self.setup_worklist(f"c{self.combination_idx}_fill_helper_plate.gwl")
        wl.add(liha.sterile_wash())
        wl.add(
            roma.move_plate(
                self.setup.overnight_12col,
                lid1_pos,
                new_lid_gridsite=plate1_pos,
                end_with_covered_plate=False,
            )
        )
        wl.add(
            roma.move_plate(
                self.setup.helper_plate,
                plate2_pos,
                new_lid_gridsite=lid2_pos,
                end_with_covered_plate=False,
            )
        )

        ## Transfer overnight cultures to 384 well helper plate using Liha
        for replicate in range(4):
            if replicate % 2:
                column_mask = self.column_mask384_B
            else:
                column_mask = self.column_mask384_A
            if replicate < 2:
                start_col = 1
            else:
                start_col = 2
            wl.add(
                liha.mix(
                    self.setup.overnight_12col,
                    self.config.overnight_culture_cols[replicate],
                    250,
                    2,
                    self.tip_arr8,
                    tip_array=self.tip_arr8,
                )
            )
            wl.add(
                liha.fill_plate(
                    self.setup.overnight_12col,
                    self.setup.helper_plate,
                    40,
                    self.column_mask96,
                    column_mask,
                    liquid_class="LB FD ZMAX",
                    tip_array=self.tip_arr,
                    skip=1,
                    src_col=self.config.overnight_culture_cols[replicate],
                    start_col=start_col,
                    end_col=24,
                ),
            )
            wl.add(liha.sterile_wash())

        wl.add(
            roma.incubate(self.setup.overnight_12col), msg="filling_helper_plate_done"
        )
        wl.add(roma.incubate(self.setup.helper_plate))

    def prefill_assay_plates(self):
        wl = self.setup_worklist(f"c{self.combination_idx}_prefill_assay.gwl")
        for assay_plate in self.setup.assay_plates:
            wl.add(
                roma.move_plate(
                    assay_plate,
                    plate2_pos,
                    new_lid_gridsite=lid2_pos,
                    end_with_covered_plate=False,
                )
            )
            # Prefill assay plates with medium
            wl.add(
                mca.fill_384plate(
                    self.setup.medium,
                    assay_plate,
                    self.config.assay_total_vol * 0.9,
                    get_tips=True,
                    liquid_class="Water free dispense",
                )
            )
            wl.add(
                roma.incubate(assay_plate),
                msg=f"assay_plate_filled",
            )

    def infect_assays(self):
        wl = self.setup_worklist(f"c{self.combination_idx}_infect_assays.gwl")
        wl.add(mca.get_pintool())
        wl.add(mca.clean_pintool())
        wl.add(mca.drop_pintool())
        wl.add(mca.get_pintool())
        wl.add(
            roma.move_plate(
                self.setup.helper_plate,
                plate2_pos,
                new_lid_gridsite=lid2_pos,
                end_with_covered_plate=False,
            )
        )
        for assay_plate in self.setup.assay_plates:
            wl.add(
                roma.move_plate(
                    assay_plate,
                    plate1_pos,
                    new_lid_gridsite=lid1_pos,
                    end_with_covered_plate=False,
                )
            )
            wl.add(
                mca.replicate_with_pintool(
                    self.setup.helper_plate, assay_plate, dip_offset=10
                )
            )
            wl.add(mca.move(self.setup.helper_plate))
            wl.add(roma.incubate(assay_plate))
        wl.add(start_timer(2))
        wl.add(mca.clean_pintool())
        wl.add(mca.drop_pintool())
        wl.add(roma.incubate(self.setup.helper_plate))

    def treat_cultures(self):
        wl = self.setup_worklist(f"c{self.combination_idx}_treat_strains.gwl")
        wl.add(wait_timer(2, self.config.exponential_growth_time))
        combination = self.combination

        # important: first handle _II then _I, as _II is lower concentrated.
        for assayplate, antibioticplate in zip(
            [combination["assay_II"], combination["assay_I"]],
            [combination["antibiotics_II"], combination["antibiotics_I"]],
        ):
            wl.add(
                roma.instert_plate_to_reader(
                    assayplate,
                    end_with_covered_plate=False,
                    new_lid_gridsite=lid2_pos,
                )
            )
            wl.add(self.setup.lumread.measure(assayplate.name + "_pre.xml"))
            wl.add(
                roma.move_plate(
                    assayplate,
                    plate2_pos,
                    end_with_covered_plate=False,
                ),
            )
            wl.add(
                roma.move_plate(
                    antibioticplate,
                    plate3_pos,
                    new_lid_gridsite=lid1_pos,
                    end_with_covered_plate=False,
                )
            )
            wl.add(mca.get_tips())
            for row, col in product([1, 2], [1, 2]):
                wl.add(
                    mca.fill_rep(
                        antibioticplate,
                        assayplate,
                        self.config.assay_total_vol * 0.1,
                        liquid_class="Water free dispense",
                        row=row,
                        col=col,
                    ),
                    msg=f"treat_r{row}_c{col}",
                )
            wl.add(mca.return_tips())
            wl.add(
                roma.instert_plate_to_reader(
                    assayplate,
                    end_with_covered_plate=False,
                    new_lid_gridsite=lid2_pos,
                )
            )
            wl.add(self.setup.lumread.measure(assayplate.name + "_i0.xml"))
            wl.add(roma.incubate(assayplate))
            wl.add(roma.store(antibioticplate))

    def luminescence_read_loop(self):
        approx_time_per_loop_min = (
            len(self.setup.assay_plates) * self.config.approx_plate_read_duration
        )
        n_loops = int(
            np.ceil(self.config.assay_duration_h * 60 / approx_time_per_loop_min)
        )

        assay_plate_names = [plate.name for plate in self.setup.assay_plates]
        counter_dict = dict(zip(assay_plate_names, [1] * len(assay_plate_names)))

        for i in range(n_loops):
            wl = self.setup_worklist(f"c{self.combination_idx}_loop_{str(i + 1)}.gwl")
            for assayplate in self.setup.assay_plates:
                wl.add(
                    roma.instert_plate_to_reader(
                        assayplate,
                        new_lid_gridsite=lid2_pos,
                        end_with_covered_plate=False,
                    )
                )
                wl.add(
                    self.setup.lumread.measure(
                        assayplate.name
                        + "_i"
                        + str(counter_dict[assayplate.name])
                        + ".xml"
                    )
                )
                counter_dict[assayplate.name] += 1
                wl.add(roma.incubate(assayplate))

    def dilute_second_overnight(self):
        # BUYS APPROX 5 HOURS TIME
        wl = self.setup_worklist(f"c{self.combination_idx}_dilute_second_overnight.gwl")
        wl.add(
            roma.move_plate(
                self.setup.next_overnight_12col,
                shaker_pos,
                new_lid_gridsite=plate3_pos,
                end_with_covered_plate=False,
            ),
            msg="prefill_lb_in_next_overnight_plate",
        )
        wl.add(liha.sterile_wash())

        wl.add(
            roma.move_plate(
                self.setup.overnight_12col,
                lid1_pos,
                new_lid_gridsite=plate1_pos,
                end_with_covered_plate=False,
            )
        )
        for i, col in enumerate(self.config.overnight_culture_cols):
            tip_arr = 8 * [False]
            tip_arr[i + 1] = True
            wl.add(
                liha.mix(
                    self.setup.overnight_12col,
                    col,
                    250,
                    2,
                    self.tip_arr8,
                    tip_array=self.tip_arr8,
                ),
            )
            wl.add(
                liha.vol_transfer(
                    self.setup.overnight_12col,
                    col,
                    self.setup.next_overnight_12col,
                    col,
                    self.config.next_on_transfer_vol,
                    tip_array=tip_arr,
                    column_mask=tip_arr,
                ),
                msg=f"transfering_strain_col{col}_to_next_overnight_plate",
            )

        wl.add(liha.sterile_wash())
        wl.add(roma.incubate(self.setup.next_overnight_12col))
        wl.add(
            roma.incubate(self.setup.overnight_12col),
            msg="dilute_second_overnight_done",
        )

    def setup_worklist(self, name):
        return self.experiment.setup_worklist(
            name, protocol=self.protocol, key="wl_II", folder_name="worklists_II"
        )
