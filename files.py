import os
import pandas as pd
import numpy as np
from itertools import product
from general_classes import PathManager


class SetupFiles:
    def __init__(self, config, experiment):
        self.experiment = experiment
        self.config = config
        self.pm = PathManager(basepath="current")
        self.folder_path = self.experiment.paths["plate_files"]
        self.drugs_df = pd.read_excel(
            self.pm.file_path("drugs.xlsx", folder="plate_files")
        )
        self.drugs = list(self.drugs_df.drug.values)
        self.ref_concentrations = dict(zip(self.drugs_df.drug, self.drugs_df.micRef))
        self.rel_concentrations = config.x
        self.reservoir_dfs = {}

    def save_csv(self, df, filename):
        path = os.path.join(self.folder_path, filename)
        df.to_csv(path, index=False)

    def save_json(self, df, filename):
        path = os.path.join(self.folder_path, filename)
        df.to_json(path, orient="records", indent=2)

    def mk_fill_table(self, V=20):
        self.drugs_df["cmax"] = self.drugs_df.micRef * max(self.config.x)
        self.drugs_df["cmix"] = self.drugs_df["cmax"] * 10 * 2
        self.drugs_df[f"V_stock per {V}ml [ml]"] = (self.drugs_df["cmix"] * V) / (
            self.drugs_df.stock * 1000
        )
        self.drugs_df[f"V_LB [ml]"] = V - self.drugs_df[f"V_stock per {V}ml [ml]"]
        self.experiment.save_csv(self.drugs_df, "fill_table.csv", folder_key="notes_I")

    def mk_reservoir_files(self):
        for antibiotic, reference_concentration in self.ref_concentrations.items():
            absolute_concentrations = (
                np.array(self.rel_concentrations)
                * reference_concentration
                * 10  # treat 1:10 and
                * 2  # mix 1:1
            )
            reservoir_df = pd.DataFrame(
                {
                    "col": list(range(1, 13)),
                    "conc": absolute_concentrations,
                }
            )
            self.save_json(reservoir_df, f"{antibiotic}_reservoir.json")
            self.reservoir_dfs[antibiotic] = reservoir_df

    @staticmethod
    def mk_rotated_reservoir(reservoir_B):
        reservoir_B_I = reservoir_B.loc[1:6, :].copy()
        reservoir_B_I["row"] = list("BCDEFG")
        reservoir_B_I.set_index("row", inplace=True)
        reservoir_B_II = reservoir_B.loc[7:12, :].copy()
        reservoir_B_II["row"] = list("BCDEFG")
        reservoir_B_II.set_index("row", inplace=True)
        reservoir_B_I.loc["A", "conc"] = None
        reservoir_B_I.loc["H", "conc"] = None
        reservoir_B_II.loc["A", "conc"] = None
        reservoir_B_II.loc["H", "conc"] = None
        return reservoir_B_I, reservoir_B_II

    @staticmethod
    def combine_drugs(reservoir_A, reservoir_B, drug_A, drug_B):
        df_I = []
        df_II = []
        reservoir_B_I, reservoir_B_II = SetupFiles.mk_rotated_reservoir(reservoir_B)
        for row, col in product("ABCDEFGH", range(1, 13)):
            well = f"{row}{col}"
            concentration_A = reservoir_A.loc[col, "conc"] / 2
            concentration_B_I = reservoir_B_I.loc[row, "conc"] / 2
            concentration_B_II = reservoir_B_II.loc[row, "conc"] / 2
            df_I.append(
                {
                    "well": well,
                    "drug_A": drug_A if pd.notna(concentration_B_I) else None,
                    "drug_B": drug_B if pd.notna(concentration_B_I) else None,
                    "conc_A": concentration_A if pd.notna(concentration_B_I) else None,
                    "conc_B": concentration_B_I,
                }
            )
            df_II.append(
                {
                    "well": well,
                    "drug_A": drug_A if pd.notna(concentration_B_I) else None,
                    "drug_B": drug_B if pd.notna(concentration_B_I) else None,
                    "conc_A": concentration_A if pd.notna(concentration_B_II) else None,
                    "conc_B": concentration_B_II,
                }
            )
        return pd.DataFrame(df_I), pd.DataFrame(df_II)

    @staticmethod
    def map_to_assay(antibiotic_plate):
        rows_96 = "ABCDEFGH"
        rows_384 = "ABCDEFGHIJKLMNOP"
        mapped_entries = []
        for _, row in antibiotic_plate.iterrows():
            source_well = row["well"]
            row_letter = source_well[0]
            col_number = int(source_well[1:])
            base_row_index = rows_96.index(row_letter)
            target_rows = [
                rows_384[2 * base_row_index],
                rows_384[2 * base_row_index + 1],
            ]
            target_cols = [2 * col_number - 1, 2 * col_number]
            for r in target_rows:
                for c in target_cols:
                    new_entry = row.copy()
                    new_entry["src_well"] = source_well
                    new_entry["well"] = f"{r}{c}"
                    mapped_entries.append(new_entry)
        df = pd.DataFrame(mapped_entries)
        df.conc_A = df.conc_A / 10
        df.conc_B = df.conc_B / 10
        return df

    def mk_combined_drug_files(self):
        for c, combination in self.config.combinations.items():
            drug_A = combination["a"]
            drug_B = combination["b"]
            reservoir_A = self.reservoir_dfs[drug_A].set_index("col")
            reservoir_B = self.reservoir_dfs[drug_B].set_index("col")
            df_I, df_II = SetupFiles.combine_drugs(
                reservoir_A, reservoir_B, drug_A, drug_B
            )
            self.save_csv(df_I, f"antibiotics_c{c}_I.csv")
            self.save_csv(df_II, f"antibiotics_c{c}_II.csv")
            assay_I = SetupFiles.map_to_assay(df_I)
            assay_II = SetupFiles.map_to_assay(df_II)
            self.save_csv(assay_I, f"assay_c{c}_I_384.csv")
            self.save_csv(assay_II, f"assay_c{c}_II_384.csv")
            combination.update(
                {
                    "Pab_I": df_I,
                    "Pab_II": df_II,
                    "assay_I": assay_I,
                    "assay_II": assay_II,
                }
            )
