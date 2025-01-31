import os
import pandas as pd
from platereads import XmlFile
import numpy as np
import sys

upper = {
    "C": True,
    "E": True,
    "G": True,
    "I": False,
    "K": False,
    "M": False,
}


class PathManager:
    def __init__(self):
        self.basepath = os.path.dirname(os.getcwd())
        self.create_folder("analysis")
        self.create_folder("analysis/obj")
        self.analysispath = self.folder_path("analysis")
        self.objpath = self.folder_path("analysis/obj")

    def folder_path(self, foldername):
        return os.path.join(self.basepath, foldername)

    def create_folder(self, foldername):
        path = self.folder_path(foldername)
        if not os.path.exists(path):
            os.makedirs(path)

    def check_files(self, names, folder):
        folder_path = self.folder_path(folder)
        existing_files = os.listdir(folder_path)
        return [any(file.startswith(name) for file in existing_files) for name in names]


class Experiment(PathManager):
    def __init__(
        self,
        treatment_keyword="antibiotic",
        load_setup=True,
        load_config=True,
        polymeasure=False,
    ):
        super().__init__()
        self.treatment_keyword = treatment_keyword
        self.polymeasure = polymeasure
        self.assay_df = pd.DataFrame()
        self.agar_df = pd.DataFrame()
        self.timelog_df = self.load_timelog()

        self.get_plate_names()

        self.load_xml_files()
        self.get_treatment_times()
        self.get_read_times()
        time_difference = self.assay_df["read_times"] - self.assay_df["treatment_times"]

        self.assay_df["treatment_duration"] = time_difference.dt.total_seconds() / 3600
        self.assay_df.dropna(inplace=True)
        self.df = self.summarize_xmls()

        if load_setup:
            self.antibiotic_p_dicts = self.load_setup()
            self.assing_setup()
            self.df["control"] = self.df.antibiotic.isnull()
            self.control = self.df[self.df.control]
            self.df = self.df[self.df.control == False]
            self.df["index"] = self.df["file_name"] + "_" + self.df.well
            self.df.set_index("index", inplace=True)
            self.antibiotic_dict = dict(zip(self.df.index, self.df.antibiotic.values))
        if load_config:
            self.load_config()

        self.time_dict = dict(
            zip(self.df.file_name.values, self.df.treatment_duration.values)
        )

    def load_config(self):
        sys.path.append(os.path.join(self.basepath, "src_code"))
        from configurations import AssayConfiguration

        self.config = AssayConfiguration()
        self.concentration_dict = {}
        for sub_dict in self.config.layout.values():
            self.concentration_dict.update(
                {sub_dict["antibiotic"]: sub_dict["working_concentration"]}
            )

    def load_setup(self):
        path = os.path.join(self.folder_path("notes"), "antibiotic_plate.xlsx")
        setup = pd.read_excel(path).dropna(how="all", axis=1)
        rows = np.array(["C", "E", "G", "I", "K", "M"])
        wells = [row + "5" for row in rows.tolist()]
        antibiotic_p_dicts = {}
        for i, col in enumerate(setup.columns[1:]):
            antibiotic_p_dicts.update({"p" + str(i): dict(zip(wells, setup[col][1:7]))})
        return antibiotic_p_dicts

    def assing_setup(self):
        for p in self.antibiotic_p_dicts:
            df = self.df[self.df.p == p]
            self.df.loc[df.index, "antibiotic"] = df.well.map(
                self.antibiotic_p_dicts[p]
            )

    def get_plate_names(self):
        self.agar_df["name"] = list(
            set(
                filter(
                    lambda x: x is not None,
                    self.timelog_df.comment.apply(
                        lambda x: "_".join(x.split("_")[1:5]) if "agar" in x else None
                    ),
                )
            )
        )

        self.assay_df["name"] = list(
            set(
                filter(
                    lambda x: (x is not "") & (x is not None),
                    self.timelog_df.comment.apply(
                        lambda x: "_".join(x.split("_")[2:5]) if "assay" in x else None
                    ),
                )
            )
        )

        self.agar_df.dropna(inplace=True)
        self.assay_df.dropna(inplace=True)
        self.assay_df["exists"] = self.check_files(self.assay_df.name, "lum_files")
        self.assay_df.reset_index(drop=True, inplace=True)
        self.agar_df.reset_index(drop=True, inplace=True)
        self.agar_df[["type", "t", "p", "d"]] = self.agar_df.name.apply(name_elements)
        self.assay_df[["type", "t", "p"]] = self.assay_df.name.apply(name_elements)

    def load_xml_files(self):
        for i, row in self.assay_df.iterrows():
            if row.exists:
                path = os.path.join(self.folder_path("lum_files"), row["name"] + ".xml")
                self.assay_df.loc[i, "xml_file"] = XmlFile(path, self.polymeasure)

    def load_timelog(self):
        filepath = os.path.join(self.folder_path("notes"), "timelog.csv")
        return pd.read_csv(filepath)

    def get_treatment_times(self):
        treatments = self.timelog_df[
            self.timelog_df.comment.str.contains(self.treatment_keyword)
        ]
        treat_dict = {}
        for i, time in enumerate(treatments.datetime):
            treat_dict.update({"p" + str(i): pd.to_datetime(time).tz_localize("CET")})
            self.assay_df["treatment_times"] = self.assay_df.p.map(treat_dict)

    def get_read_times(self):
        times = {}
        for i, row in self.assay_df.iterrows():
            if row.exists:
                times.update(
                    {
                        row["name"]: pd.to_datetime(
                            row.xml_file.end_datetimes[0]
                        ).tz_convert("CET")
                    }
                )
        self.read_times = times
        self.assay_df["read_times"] = self.assay_df.name.map(times)

    def summarize_xmls(self):
        df = []
        for _, row in self.assay_df.iterrows():
            sub_df = row.xml_file.df[
                ["well", "row", "column", "signal", "file_name"]
            ].copy()
            sub_df.loc[
                :, ["treatment_times", "read_times", "treatment_duration", "t", "p"]
            ] = (
                row.treatment_times,
                row.read_times,
                row.treatment_duration,
                row.t,
                row.p,
            )
            df.append(sub_df)
        return pd.concat(df).reset_index()

    def map_antibiotics(self, layout):
        self.df["antibiotic"] = self.df.row.map(upper)
        self.df.antibiotic = self.df.apply(
            lambda row: layout[row.p][row.antibiotic], axis=1
        )


def disect_name(name):
    elements = name.split("_")
    row = {"name": name, "type": elements[0], "t": elements[1], "p": elements[2]}
    if len(elements) > 3:
        row.update({"d": elements[3]})
    return pd.Series([row])


def name_elements(name):
    elements = name.split("_")
    row = {"type": elements[0], "t": elements[1], "p": elements[2]}
    if len(elements) > 3:
        row["d"] = elements[3]
    return pd.Series(row)
