from dataclasses import dataclass, field
from general_classes import PathManager
import pandas as pd
from typing import Dict, List, Optional, Tuple
import itertools

pm = PathManager(basepath="current")


drugs_df = pd.read_excel(pm.file_path("drugs.xlsx", folder="plate_files"))


@dataclass(frozen=True)
class AssayConfiguration:
    # Pi - drug reservoir gradient plate
    x: List[float] = field(
        default_factory=lambda: [64, 32, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0]
    )
    columns: List[int] = field(default_factory=lambda: list(range(1, 13)))
    antibiotic_reservoir_plate_final_vol: int = 10000
    pi_feeding_cols: Dict[int, Optional[int]] = field(init=False)
    pi_column_concentrations: Dict[int, float] = field(init=False)
    pi_transfer_volumes: Dict[int, Optional[float]] = field(init=False)
    pi_LB_fill_volumes: Dict[int, float] = field(init=False)  # New dictionary
    drugs: List[str] = field(init=False)  # Now an instance variable

    # important to not take too much. 6x100 = 600 μL (which is fine if 1 mL is in the source well.)
    # 100 μL is fine as we only need 4×6 for the assay plate
    drug_plate_vol: int = 100

    # subreservoir
    drug_sub_reservoir_vol: int = 125

    # assayplate
    assay_total_vol: int = 60

    # durations
    # use 10ml for each rep
    overnight_incubation_time: int = 12 * 60 * 60  # 14 hours
    exponential_growth_time: int = int(2 * 60 * 60)  # 2 hours
    assay_duration_h: int = 5  # 5 hours
    approx_plate_read_duration: int = 5  # minutes

    # overnight plate (12 col)
    strain_reservoir_well_volume: int = 10000
    overnight_culture_cols: List[int] = field(init=False)
    combinations: Dict[Tuple[str, str], Dict[str, Optional[float]]] = field(init=False)

    # next overnight plate (12 col)
    next_on_vol = 5000
    next_on_transfer_vol = 5

    def __post_init__(self):
        object.__setattr__(
            self,
            "pi_feeding_cols",
            dict(zip(self.columns, [None, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, None])),
        )
        object.__setattr__(
            self, "pi_column_concentrations", dict(zip(self.columns, self.x * 2))
        )
        object.__setattr__(
            self, "pi_transfer_volumes", self._calculate_pi_transfer_volumes()
        )
        object.__setattr__(
            self, "pi_LB_fill_volumes", self._calculate_pi_LB_fill_volumes()
        )  # New fill volume calculation
        object.__setattr__(
            self, "drugs", list(drugs_df.drug.values)
        )  # Ensure drugs is an instance variable
        object.__setattr__(self, "overnight_culture_cols", [1, 4, 7, 10])
        index_pairs = list(itertools.combinations(range(len(self.drugs)), 2))
        combo_dict = {}
        counter = 0
        for idx1, idx2 in index_pairs:
            drug1 = self.drugs[idx1]
            drug2 = self.drugs[idx2]
            combo_dict[counter] = {
                "i": idx1,
                "j": idx2,
                "a": drug1,
                "b": drug2,
            }
            counter += 1
        object.__setattr__(self, "combinations", combo_dict)

    def _calculate_pi_transfer_volumes(self) -> Dict[int, Optional[float]]:
        transfer_volumes = {}
        for target_col, src_col in self.pi_feeding_cols.items():
            if src_col is None:
                transfer_volumes[target_col] = None
                continue
            x_target = self.pi_column_concentrations[target_col]
            x_src = self.pi_column_concentrations[src_col]
            transfer_vol = self.antibiotic_reservoir_plate_final_vol * x_target / x_src
            transfer_volumes[target_col] = transfer_vol
        return transfer_volumes

    def _calculate_pi_LB_fill_volumes(self) -> Dict[int, float]:
        fill_volumes = {}
        for col in self.columns:
            transfer_vol = self.pi_transfer_volumes.get(col, 0) or 0  # None counts as 0
            fill_volumes[col] = self.antibiotic_reservoir_plate_final_vol - transfer_vol

        # Explicitly set fill volume for column 1 to 0
        fill_volumes[1] = 0
        return fill_volumes

    def __repr__(self) -> str:
        attributes = {
            attr: getattr(self, attr)
            for attr in dir(self)
            if not attr.startswith("__") and not callable(getattr(self, attr))
        }
        formatted_attrs = ",\n    ".join(f"{k} = {v}" for k, v in attributes.items())
        return f"AssayConfiguration(\n    {formatted_attrs}\n)"
