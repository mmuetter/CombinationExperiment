import pandas as pd
import os
import numpy as np


class PlateFiles:
    def __init__(self, config, experiment):
        self.path = experiment.paths["notes"]
        self.antibiotics, self.concentrations = self.create_antibiotic_plate()
        self.save_sheets_to_excel(
            [self.antibiotics, self.concentrations],
            ["antibiotic", "concentration"],
            "antibiotic_plate.xlsx",
        )
        self.fill_table = self.create_fill_table()
        self.save_sheets_to_excel(
            [self.fill_table], ["mix_instructions"], "fill_table.xlsx"
        )

    def create_antibiotic_plate(self):
        antibiotics = pd.DataFrame(
            np.full((8, 12), None), index=list("ABCDEFGH"), columns=range(1, 13)
        )
        concentrations = pd.DataFrame(
            np.full((8, 12), None), index=list("ABCDEFGH"), columns=range(1, 13)
        )
        row_dict = {"a": ["B", "C", "D"], "b": ["E", "F", "G"]}
        col_dict = {"I": 1, "II": 3}
        for key, info in self.layout.items():
            column = col_dict[key[:-1]]
            rows = row_dict[key[-1]]
            for row in rows:
                antibiotics.loc[row, column] = info["antibiotic"]
                concentrations.loc[row, column] = info["working_concentration"] * 10
        return antibiotics, concentrations

    def create_fill_table(self):
        df = []
        for value in self.layout.values():
            df.append(value)
        df = pd.DataFrame().from_records(df)
        df["10xc0"] = df.working_concentration * 10
        df["Vs[ul] per 1ml medium"] = df["10xc0"] / df.stock_concentration
        return df

    def save_sheets_to_excel(self, sheets: list, sheet_names: list, name):
        with pd.ExcelWriter(os.path.join(self.path, name)) as writer:
            for sheet_name, sheet in zip(sheet_names, sheets):
                sheet.to_excel(writer, sheet_name=sheet_name)
