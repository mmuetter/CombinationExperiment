import pandas as pd
from figures import Figure, styles
from general_classes import PathManager
from pypetting_experiments.luminescence_validation.eval_classes import plot_df

pm = PathManager()

figure = Figure(style=styles["paper"])
figure.add_folder(pm.folder_path("analysis"))

df = pd.read_csv("df.csv")
antibiotics = df.antibiotic.unique()

for antibiotic in antibiotics:
    _, ax = figure.create_figure_with_style()
    plot_df(df, antibiotic, ax1=ax)
    figure.save_figure(antibiotic + "_comp")
