import pandas as pd
from general_classes import PathManager
from pypetting_experiments.luminescence_validation.eval_classes import plot_df
from luminescence_paper_settings import figure_si_double as si_figure
from luminescence_paper_settings import figure_main_single as main_figure
from luminescence_paper_settings import color_dict, markers_dict

pm = PathManager()


df = pd.read_csv("df.csv")
antibiotics = df.antibiotic.unique()

main = []
si = ["Mecilinam", "Piperacillin"]


###########################

sub_df = df[df.antibiotic == antibiotic]


df_pivot = sub_df.pivot_table(
    index=[col for col in sub_df.columns if col not in ["signal_type", "signal"]],
    columns="signal_type",
    values="signal",
    aggfunc="first",  # Handle duplicates if any exist
).reset_index()

df_merged = df_pivot.groupby(["ti", "well"], as_index=False).first()


df_pivot.columns = [
    col[0] if isinstance(col, tuple) else col for col in df_pivot.columns
]

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.pointplot(data=df_merged, x="cfu", y="rlu", hue="ti")
g.set(yscale="log", xscale="log")
plt.show()


###########################

for antibiotic in antibiotics:
    if antibiotic in main:
        figure = main_figure
        _, ax = figure.create_figure_with_style()
        plot_df(
            df,
            antibiotic,
            ax=ax,
            lum_signal_type="rlu",
            lum="signal",
            s=25,
            rlu_color=color_dict["lum"],
            cfu_color=color_dict["cfu"],
            cfu_marker=markers_dict["cfu"],
            rlu_marker=markers_dict["lum"],
        )
        figure.save_figure(antibiotic + "_comp")
        figure.save_figure(antibiotic + "_comp", folderpath=pm.folder_path("analysis"))

    if antibiotic in si:
        figure = si_figure
        _, ax = figure.create_figure_with_style()
        plot_df(
            df,
            antibiotic,
            ax=ax,
            lum_signal_type="rlu",
            lum="signal",
            s=25,
            rlu_color=color_dict["lum"],
            cfu_color=color_dict["cfu"],
            cfu_marker=markers_dict["cfu"],
            rlu_marker=markers_dict["lum"],
        )
        figure.save_figure(antibiotic + "_comp")
        figure.save_figure(antibiotic + "_comp", folderpath=pm.folder_path("analysis"))
