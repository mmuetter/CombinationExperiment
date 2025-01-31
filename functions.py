import pandas as pd
import numpy as np


def layout_antibiotic_plate(drugs, pA, pB):
    if pA + pB != 1:
        raise ValueError("The sum of the proportions must be 1.")

    antibiotics = drugs.drug
    rows = drugs.row
    cols = list(np.array(range(len(drugs))) + 2)
    row_ab = dict(zip(rows, antibiotics))
    col_ab = dict(zip(cols, antibiotics))

    well_list = []
    for row in row_ab.keys():
        for col in col_ab.keys():
            well = f"{row}{col}"
            well_list.append([well, row, col, row_ab[row], col_ab[col]])

    df = pd.DataFrame(
        well_list, columns=["well", "row", "col", "antibiotic_A", "antibiotic_B"]
    )
    df.set_index("well", inplace=True)
    df["pA"] = pA
    df["pB"] = pB
    return df
