from eval_agar_plates import SectionManager, ImManager
import os
from general_classes import TimeLog
import numpy as np
import pandas as pd


def recover_sec_df(path_manager):
    timelog = TimeLog(os.path.dirname(os.getcwd()))
    timelog.set_t0("antibiotic_added_to_deepwell")
    img_manager = ImManager("agar_t0_p0_d4_n0.png", timelog)
    key = "spot_agar"
    img_manager.spot_times = img_manager.timelog.get_timepoints(key)
    keys = [
        img_manager.extract_name(key, 2, 5)
        for key in list(img_manager.spot_times.keys())
    ]
    img_manager.spot_times_h = dict(
        zip(keys, np.array(list(img_manager.spot_times.values())) / 3600)
    )
    identifier = set(
        [img_manager.identifier(key) for key in list(img_manager.spot_times.keys())]
    )
    section_manager = SectionManager(identifier, img_manager.path_manager)
    section_manager.sec_df = pd.read_csv(
        path_manager.file_path("sections.csv", "analysis/obj"), index_col=False
    )

    section_manager.sec_df["picked_sec_name"] = section_manager.sec_df["sec_name"]
    section_manager.sec_df["sec"] = section_manager.sec_df.sec_id.apply(
        lambda x: x.split("_")[-1]
    )
    section_manager.sec_df["prefix"] = section_manager.sec_df.sec_id.apply(
        lambda x: "_".join(x.split("_")[:-1])
    )
    section_manager.sec_df["dilution"] = section_manager.sec_df.picked_sec_name.apply(
        lambda x: x.split("_")[3]
    )
    section_manager.sec_df.to_csv(section_manager.filepath)
