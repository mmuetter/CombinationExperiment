from platereads import Plate, Setup, LumData
import os
from general_classes import Setup


class ImportLumMeasures:
    def __init__(self, path=os.path.dirname(os.getcwd())):
        self.path = path
        self.assay_row_dict384 = {
            chr(i): (i - 64) // 2 for i in range(ord("C"), ord("N"), 2)
        }
        self.assay_row_dict96 = {chr(i): i - 65 for i in range(ord("B"), ord("H"))}
        self.load_plate()

    def load_plate(self):
        self.plate = Plate(
            "",
            path=self.path,
            filetype=".xml",
            polymeasure=False,
        )
        self.raw_data = self.plate.file_summary.copy()

    def add_setup(self, setup: Setup, well_col="well"):
        setup.map_wells(self.raw_data, well_col=well_col)

    def set_tstart(self, timelog):
        self.raw_data["time"] = self.raw_data.file_time_start - timelog.t_start
        self.t_start = timelog.t_start
        self.raw_data["t"] = self.raw_data.time.dt.total_seconds() / 3600

    def remove_background(self, indicator_col):
        self.raw_data.loc[:, "control"] = self.raw_data[indicator_col].isnull()
        lumdata = LumData(self.raw_data)
        lumdata.closest_control_noise()
        self.df = lumdata.signal_data
