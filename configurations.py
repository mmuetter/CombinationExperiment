from dataclasses import dataclass, field
import numpy as np
from pypetting.labware import labwares
from pypetting import Labware


@dataclass(frozen=True)
class AssayConfiguration:
    concentration_gradient = [0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32]
    plate_names = [
        "antibiotics32",
        "antibiotics16",
        "antibiotics8",
        "antibiotics4",
        "antibiotics2",
        "antibiotics1",
        "antibiotics05",
        "antibiotics025",
        "antibiotics0125",
        "antibiotics0",
    ]

    antibiotic_reservoir_plate_final_vol = 1000

    # important to not take to much. 6x100 = 600 ul (which is fine if 1ml is in the src well.)
    # 100 ul is fine as we only need 4x10 for the assay plate
    drug_plate_vol = 100

    ## assayplate
    assay_total_vol = 60
    assay_col_start = 3
    assay_col_end = 15

    # durations
    overnight_incubation_time = 16 * 60 * 60  # 16 hours
    exponential_growth_time = 1.5 * 60 * 60  # 90 minutes
    assay_duration = 4 * 60 * 60  # 4 hours
    approx_plate_read_duration = 3  # minutes

    # overnight plate (12 col)
    # should last for 4 plates. 5 * 6 * 6 * 50 ul = 9000 ul
    strain_reservoir_well_volume = 10000
    overnight_culture_cols = [1, 4, 7, 10]
    exponential_culture_cols = [2, 5, 8, 11]
