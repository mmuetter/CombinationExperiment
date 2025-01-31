import os
from eval_agar_plates import (
    ImManager,
)
from general_classes import TimeLog, Setup
import pandas as pd
from pypetting_experiments.luminescence_validation.eval_classes import (
    ImportLumMeasures,
)


def full_evaluation(path_manager, background_img_name="agar_t4_p0_d3_n0.png"):
    path_manager.make_dir("analysis")
    path = os.path.dirname(os.getcwd())

    timelog = TimeLog(path)
    timelog.set_t0("antibiotic_added_to_deepwell")

    setup = Setup(path_manager.file_path("antibiotic_plate.xlsx", "notes"))

    src_well_map = {
        "p0": {
            "C5": "B1",
            "E5": "C1",
            "G5": "D1",
            "I5": "E1",
            "K5": "F1",
            "M5": "G1",
        },
        "p1": {
            "C5": "B3",
            "E5": "C3",
            "G5": "D3",
            "I5": "E3",
            "K5": "F3",
            "M5": "G3",
        },
    }

    def col_map_function(self):
        self.sec_df["plate"] = self.sec_df.prefix.apply(
            lambda x: self.extract_name(x, 2, 3)
        )
        self.sec_df["col"] = 5
        self.sec_df["well"] = self.sec_df.apply(lambda x: x.row + str(x.col), axis=1)
        for plate in self.sec_df.plate.unique():
            print(src_well_map[plate])
            mask = self.sec_df.plate == plate
            self.sec_df.loc[mask, "antibiotic_src_well"] = self.sec_df.loc[
                mask, "well"
            ].map(src_well_map[plate])

    rec_settings = {
        "max_pixels": 2000,
        "disk_size": 3,
        "low_threshold": 20,
        "distance_coefficient": 0.6,
        "min_pixels": 15,
        "overlap_threshold": 0.75,
        "min_radius": 3,
    }

    img_manager = ImManager(background_img_name, timelog, absorbance_const=0.4)
    img_manager.process_images()
    img_manager.select_sections(redo=False)
    img_manager.map_wells(col_map_function)
    img_manager.evaluate_colonies(redo=False, rec_settings=rec_settings)
    img_manager.add_setup(setup)
    img_manager.map_setup(well_col="antibiotic_src_well")

    # img_manager.reeval_section("t3_p0_sec2", rec_settings=rec_settings)

    lum = ImportLumMeasures(path_manager.folder_path("lum_files"))

    lum.raw_data["plate"] = lum.raw_data.file_name.apply(lambda x: x.split("_")[-1])

    for plate in lum.raw_data.plate.unique():
        mask = lum.raw_data.plate == plate
        lum.raw_data.loc[mask, "antibiotic_src_well"] = lum.raw_data.loc[
            mask, "well"
        ].map(src_well_map[plate])

    lum.add_setup(setup, well_col="antibiotic_src_well")
    lum.set_tstart(timelog)
    lum.remove_background("antibiotic")
    lum.df.rename(columns={"method": "signal_type"}, inplace=True)

    df = pd.concat([lum.df, img_manager.df])
    df.to_csv("df.csv")
    return img_manager, lum
