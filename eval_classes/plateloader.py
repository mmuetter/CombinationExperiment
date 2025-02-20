import re
import numpy as np
from general_classes import PathManager
import pandas as pd
import os
from platereads import Plate


class PlateLoader:
    def __init__(
        self, pm: PathManager, plate_name: str, experiment: str, prep_folder: str
    ) -> None:
        self.pm = pm
        self.plate_name = plate_name
        self.experiment = experiment
        self.prep_folder = prep_folder
        self.rel_conc = self._extract_number(plate_name)
        self.exp_pm = PathManager(self.pm.folder_path(experiment))
        self.plate = self.load_plate()
        self.setup_df96 = self.load_setup()
        self.setup_df = self._map96_to384(self.setup_df96)
        self.drugs = list(
            set(
                self.setup_df.antibiotic_A.unique().tolist()
                + self.setup_df.antibiotic_B.unique().tolist(),
            )
        )
        self._calculate_rel_conc(self.setup_df)
        features = self._mk_features()
        self.df = self._map_setup(self.plate.file_summary, self.setup_df, features)
        self.df["treatment_duration"] = (
            self.df["file_time_start"] - self.df["file_time_start"].min()
        )
        self.df["t"] = self.df["treatment_duration"].dt.total_seconds() / 3600

    def _map96_to384(self, df: pd.DataFrame) -> pd.DataFrame:
        new_rows = []
        for _, row in df.iterrows():
            # Convert row letter (A..H) to integer index 0..7
            i = ord(row["row"]) - ord("A")
            # For columns 1..12
            c = int(row["col"])

            # Expanded rows will be row_letters[0], row_letters[1]
            row_letters = [chr(ord("A") + 2 * i), chr(ord("A") + 2 * i + 1)]
            # Expanded columns will be 2*c - 1, 2*c
            cols = [2 * c - 1, 2 * c]

            for r in row_letters:
                for col in cols:
                    new_entry = row.copy()
                    new_entry["row"] = r
                    new_entry["col"] = col
                    new_entry["wells"] = f"{r}{col}"
                    new_rows.append(new_entry)
        df384 = pd.DataFrame(new_rows)
        return self._assign_replicates(df384)

    def _assign_replicates(self, df: pd.DataFrame) -> pd.DataFrame:
        # Convert row letter (A..P) to an index (0..15)
        row_idx = df["row"].apply(lambda x: ord(x) - ord("A"))
        # Convert col number to zero-based index
        col_idx = df["col"].astype(int) - 1
        # Top-left = replicate 1, bottom-left = 2, top-right = 3, bottom-right = 4
        df["replicate"] = 2 * (col_idx % 2) + (row_idx % 2) + 1
        return df

    def _mk_features(
        self,
        features=[
            "replicate",
            "antibiotic_A",
            "antibiotic_B",
            "pA",
            "pB",
            "rel_conc",
            "combination_id",
        ],
    ):
        for antibiotic in self.drugs:
            features.append("f_" + antibiotic)
            features.append("rel_conc_" + antibiotic)
        return features

    def _calculate_rel_conc(self, setup_df: pd.DataFrame) -> pd.DataFrame:
        for _, row in setup_df.iterrows():
            self._write_rel_conc(setup_df, row)

    def _write_rel_conc(self, setup_df, row):
        fA = row.pA * (np.array(self.drugs) == row["antibiotic_A"])
        fB = row.pB * (np.array(self.drugs) == row["antibiotic_B"])
        f = fA + fB
        for antibiotic, fi in zip(self.drugs, f):
            setup_df.loc[row.name, "f_" + antibiotic] = fi
            setup_df.loc[row.name, "rel_conc_" + antibiotic] = fi * row.rel_conc
        self.setup_df.loc[row.name, "combination_id"] = self._name_combination(f)

    def _name_combination(self, f):
        F = np.array(f)[f.astype(bool)].tolist()
        drugs = np.array(self.drugs)[f.astype(bool)].tolist()
        name = []
        for fi, drug in zip(F, drugs):
            if fi == 1:
                fstr = "1"
            else:
                fstr = f"{fi}".replace(".", "")
            name.append(f"{drug}{fstr}")
        return "_".join(name)

    def _extract_number(self, plate_name: str) -> str:
        m = re.search(r"\d+", plate_name)
        return m.group() if m else None

    def load_plate(self) -> Plate:
        return Plate(
            self.plate_name, filetype=".xml", path=self.exp_pm.file_path("lum_files")
        )

    def load_setup(self) -> pd.DataFrame:
        setup_path = self.pm.file_path(
            f"setup_antibiotics{self.rel_conc}.csv",
            folder=os.path.join(self.prep_folder, "notes"),
        )
        df = pd.read_csv(setup_path)
        df["row"] = df["row"].str.upper()
        df["wells"] = df["row"] + df["col"].astype(str)
        return df

    def _mk_combi_id(self, setup_df):
        pA_tmp = setup_df.pA.astype(str).str.replace(".", "")
        pB_tmp = setup_df.pB.astype(str).str.replace(".", "")
        setup_df["combi_id"] = (
            setup_df["antibiotic_A"] + pA_tmp + "_" + setup_df["antibiotic_B"] + pB_tmp
        )
        return setup_df

    def _map_setup(
        self,
        df: pd.DataFrame,
        setup_df: pd.DataFrame,
        features=["antibiotic_A", "antibiotic_B", "pA", "pB", "rel_conc"],
    ) -> pd.DataFrame:
        for feat in features:
            df[feat] = df["well"].map(setup_df.set_index("wells")[feat])
        return df

    def save_df(self):
        self.pm.make_dir("dataframes")
        self.df.to_csv(
            self.pm.file_path(f"{self.plate_name}df.csv", folder="dataframes"),
            index=False,
        )
