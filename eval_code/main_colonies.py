import matplotlib.pyplot as plt
import os
import seaborn as sns
from eval_agar_plates import (
    ImManager,
)
from general_classes import TimeLog, PathManager, Setup
from figures import Figure, styles, generate_style
import seaborn as sns
import pandas as pd
from pypetting_experiments.luminescence_validation.eval_classes import (
    ImportLumMeasures,
    plot,
)


path_manager = PathManager()
path_manager.make_dir("analysis")
path = os.path.dirname(os.getcwd())

timelog = TimeLog(path)
timelog.set_t0("antibiotic_added_to_deepwell")

setup = Setup(path_manager.file_path("antibiotic_plate.xlsx", "notes"))


background_img_name = "agar_t5_p0_d1_n1.png"

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

# img_manager.reeval_section("t1_p0_sec3", rec_settings=rec_settings)

lum = ImportLumMeasures(path_manager.folder_path("lum_files"))

lum.raw_data["plate"] = lum.raw_data.file_name.apply(lambda x: x.split("_")[-1])

for plate in lum.raw_data.plate.unique():
    mask = lum.raw_data.plate == plate
    lum.raw_data.loc[mask, "antibiotic_src_well"] = lum.raw_data.loc[mask, "well"].map(
        src_well_map[plate]
    )
lum.set_tstart(timelog)
lum.add_setup(setup, well_col="antibiotic_src_well")
lum.remove_background("antibiotic")
lum.df.rename(columns={"method": "signal_type"}, inplace=True)


_, axs = plt.subplots(1, 2, figsize=(12, 16))
# fig, ax = paperfig.create_figure_with_style()
plot(lum, img_manager, "Tetracycline", axs[0])
axs[0].set_title("Tetracycline")


# fig, ax = paperfig.create_figure_with_style()
plot(lum, img_manager, "Trimethoprim", axs[1])
axs[1].set_title("Trimethoprim")

plt.show()
