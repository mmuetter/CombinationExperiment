from setup_pdcurve import (
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
)
from dataclasses import dataclass
import numpy as np
from pypetting import start_timer, wait_timer
from icecream import ic


class pdWorkflow:
    def __init__(self, setup: dataclass):
        self.experiment = experiment = setup.experiment
        self.setup = setup
        self.config = setup.config
        self.protocol = experiment.setup_protocol()
        self.tip_arr = [False] + 6 * [True] + [False]
        self.tip_arr8 = 8 * [True]
        self.column_mask96_8 = 8 * [True]
        self.column_mask96 = self.tip_arr
        self.column_mask384_A = np.array(2 * [False] + 6 * [True, False] + 2 * [False])
        self.column_mask384_B = np.array(2 * [False] + 6 * [False, True] + 2 * [False])

    def grow_overnight(self, wl_name="grow_overnight.gwl"):
        wl = self.setup_worklist(wl_name)
        wl.add(start_timer(1))
        wl.add(wait_timer(1, self.config.overnight_incubation_time))

    def prepare_exp_cultures(self):
        wl = self.setup_worklist("prepare_exponential_cultures.gwl")
        tip_vol_medium = (
            self.config.strain_reservoir_well_volume * 0.99 / sum(self.tip_arr8)
        )
        tip_transfer_vol = (
            self.config.strain_reservoir_well_volume * 0.01 / sum(self.tip_arr8)
        )

        wl.add(liha.sterile_wash())
        for exp_col in self.config.exponential_culture_cols:
            wl.add(
                liha.vol_transfer(
                    self.setup.medium_trough,
                    1,
                    self.setup.overnight_12col,
                    exp_col,
                    tip_vol_medium,
                    column_mask=self.column_mask96_8,
                    tip_array=self.tip_arr8,
                    liquid_class="Minimal CD ZMAX FAST",
                )
            )

        # Dilute Strains 1:100
        for replicate, (replicate_col, exponential_col) in enumerate(
            zip(
                self.config.overnight_culture_cols, self.config.exponential_culture_cols
            )
        ):
            wl.add(
                liha.mix(
                    self.setup.overnight_12col,
                    replicate_col,
                    250,
                    2,
                    self.column_mask96_8,
                    tip_array=self.tip_arr8,
                )
            )
            wl.add(
                liha.vol_transfer(
                    self.setup.overnight_12col,
                    replicate_col,
                    self.setup.overnight_12col,
                    exponential_col,
                    tip_transfer_vol,
                    column_mask=self.column_mask96_8,
                    tip_array=self.tip_arr8,
                    liquid_class="Minimal CD ZMAX FAST",
                )
            )
            wl.add(
                liha.mix(
                    self.setup.overnight_12col,
                    replicate_col + 1,
                    250,
                    2,
                    self.column_mask96,
                    tip_array=self.tip_arr,
                )
            )
            self._fill_replicate(replicate, wl)
            wl.add(liha.simple_wash())

        for plate in self.setup.assayplates:
            wl.add(roma.incubate(plate), msg="incubate_exp_cultures_" + plate.name)
        start_timer(2)
        wait_timer(2, self.config.exponential_growth_time)

    def _fill_replicate(self, replicate, wl):
        wl = self.setup_worklist("fill_replicate_" + str(replicate) + ".gwl")
        if replicate % 2:
            column_mask = self.column_mask384_B
        else:
            column_mask = self.column_mask384_A

        for plate in self.setup.assayplates:
            wl.add(
                roma.move_plate(
                    plate,
                    plate2_pos,
                    new_lid_gridsite=lid2_pos,
                    end_with_covered_plate=False,
                ),
                msg="fill_replicate_" + str(replicate) + "_" + plate.name,
            )
            wl.add(
                liha.fill_plate(
                    self.setup.overnight_12col,
                    plate,
                    self.config.assay_total_vol * 0.9,
                    self.column_mask96,
                    column_mask,
                    liquid_class="LB CD ZMAX FAST",
                    tip_array=self.tip_arr,
                    skip=1,
                    start_col=self.config.assay_col_start,
                    end_col=self.config.assay_col_end,
                )
            )

            wl.add(roma.store(plate))

    def treat_cultures(self):
        wl = self.setup_worklist("treat_strains.gwl")
        for assayplate, antibioticplate in zip(
            self.setup.assayplates, self.setup.antibiotic_plates
        ):
            wl.add(
                roma.move_plate(
                    assayplate,
                    plate2_pos,
                    new_lid_gridsite=lid2_pos,
                    end_with_covered_plate=False,
                ),
            )
            wl.add(
                roma.move_plate(
                    antibioticplate,
                    plate1_pos,
                    new_lid_gridsite=lid1_pos,
                    end_with_covered_plate=False,
                )
            )
            wl.add(
                mca.fill_384plate(
                    antibioticplate,
                    assayplate,
                    self.config.assay_total_vol * 0.1,
                    get_tips=True,
                ),
                msg="treat_" + assayplate.name,
            )
            wl.add(
                roma.instert_plate_to_reader(assayplate, end_with_covered_plate=False)
            )
            wl.add(self.setup.lumread.measure(assayplate.name + "_i0.xml"))
            wl.add(roma.incubate(assayplate))
            wl.add(roma.store(antibioticplate))

    def luminescence_read_loop(self, read_order, approx_time_per_read=5 * 60):
        approx_time_per_loop = len(read_order) * approx_time_per_read
        n_loops = int(np.ceil(self.config.assay_duration // approx_time_per_loop))

        counter_dict = dict(
            zip(self.setup.assay_plate_names, [0] * len(self.setup.assay_plate_names))
        )
        for i in range(n_loops):
            wl = self.setup_worklist("loop_" + str(i + 1) + ".gwl")
            for plate in read_order:
                print(counter_dict[plate])
                counter_dict[plate] += 1
                assayplate = self.setup.assay_dict[plate]
                wl.add(
                    roma.instert_plate_to_reader(
                        assayplate,
                        new_lid_gridsite=lid2_pos,
                        end_with_covered_plate=False,
                    )
                )
                wl.add(
                    self.setup.lumread.measure(
                        plate + "_i" + str(counter_dict[plate]) + ".xml"
                    )
                )
                wl.add(roma.incubate(assayplate))

    def setup_worklist(self, name):
        return self.experiment.setup_worklist(name, protocol=self.protocol)
