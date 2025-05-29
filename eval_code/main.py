import matplotlib.pyplot as plt
import os
from general_classes import TimeLog, PathManager
import os
import sys
import numpy as np


main_path = "/Users/malte/ETH-Documents/"
if main_path not in sys.path:
    sys.path.insert(0, main_path)
from combination_project import (
    figure_main_draft,
    ImportCombinationData,
    FitRates,
    PdFitter,
    plot_plate,
)

blank_path = (
    "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/noise/blank_map.json"
)
dist_ratio_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/noise/dist_ratio_map.json"

path_manager = PathManager()
path_manager.make_dir("analysis")
path = os.path.dirname(os.getcwd())
combinations = path_manager.import_json("combinations.json", "notes_III")


keys = np.array(list(combinations.keys()))
timelog = TimeLog(path, folder="notes_III")
redo_data_import = False
redo_fitting = False
summaries = {}
combination_objs = {}

for i in [5, 7, 8, 13]:  # 5, 7, 13
    key = keys[i]
    print(key)
    comb_dict = combinations[key]
    if redo_data_import:
        combination = ImportCombinationData(
            i, comb_dict, path_manager, timelog, dist_ratio_path
        )
        path_manager.save_csv(combination.plate_df, f"{combination.combi}.csv", "obj")

    comb_obj = FitRates(comb_dict, path_manager)
    if redo_fitting:
        print(f"redo fitting {key}")
        summary = comb_obj.eval_all_wells()
    else:
        summary = comb_obj.load_summary()
    summaries[key] = summary
    combination_objs[key] = comb_obj


obj = combination_objs["chloramphenicol_tetracycline"]
idx_pair = (0, 0)
fig, ax = obj.plot_concentration_idx_pair(idx_pair)
obj.plot_fit_idx_pair(idx_pair, ax=ax)
plt.show()
obj.median_heatmap(save_local=True)


obj = combination_objs["chloramphenicol_polymyxinB"]
idx_pair = (0, 0)
fig, ax = obj.plot_concentration_idx_pair(idx_pair)
obj.plot_fit_idx_pair(idx_pair, ax=ax)
plt.show()
obj.median_heatmap(save_local=True)

obj.plot_plate("I")
obj.plot_plate("II")
obj.plot_well_id("D19_I")


idx_pair = (0, 9)
obj = combination_objs["chloramphenicol_colistin"]
fig, ax = obj.plot_concentration_idx_pair(idx_pair)
obj.plot_fit_idx_pair(idx_pair, ax=ax)
plt.show()
obj.median_heatmap(save_local=True)

obj.plot_plate("I")
obj.plot_plate("II")
obj.plot_well_id("D19_I")


idx_pair = (11, 9)
obj = combination_objs["fosfomycin_tetracycline"]
fig, ax = obj.plot_concentration_idx_pair(idx_pair)
obj.plot_fit_idx_pair(idx_pair, ax=ax)
plt.show()
obj.median_heatmap(save_local=True)


df_pair = obj.idx_df(idx_pair)
display(
    df_pair[
        [
            "well",
            "t",
            "lum",
            "drop_magnitude",
            "exceed_threshold",
            "exclude_after_drop",
            "include",
        ]
    ]
)


df_pair
obj.concentrations_A
obj.concentrations_B[11]


summary
import seaborn as sns

summary_median = (
    obj.summary[["conc_A", "conc_B", "psi"]]
    .groupby(["conc_A", "conc_B"])
    .median()
    .reset_index()
)

heatmap_data = summary_median.pivot(index="conc_B", columns="conc_A", values="psi")
heatmap_data = heatmap_data.sort_index(ascending=False)  # highest conc_B at top

ax = sns.heatmap(data=heatmap_data, cmap="viridis")
ax.set_xlabel("conc_A")
ax.set_ylabel("conc_B")
plt.show()


obj.plot_heatmap()
obj.plot_topography(sigma=5)
obj.plot_smoothed_heatmap()


#######################################################
######### PD CURVES
angle_colors = {np.pi / 4: "blue", 0: "red", np.pi / 2: "green"}

from matplotlib.lines import Line2D

summary = summaries["fosfomycin_tetracycline"]
fit = PdFitter(summary)
fit.fit_psi_max()
fit.fit_mono("A")
fit.fit_mono("B")
fit.report()

phi = np.pi / 4
fit.loewe(phi, 1)
fit.bliss(phi, 30)
fit.estimate_combi_tau(phi, 1)

fit.plot_combi_pd(phi=phi, loewe=True, bliss=True)

fit.plot_combi_pd(phi=0, loewe=True, bliss=True)

fit.plot_combi_pd(phi=np.pi / 2, loewe=True, bliss=True)


fit = FitPDCurve(summary, path_manager)
_, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("r [\mu g/ml]")
ax.set_ylabel(r"$\psi [h^-1]$")
labels = []
for _ in range(200):
    df = fit.bootstrap_summary()

    # For each curve, compute radial distance in original scale and plot the points
    for phi, color in angle_colors.items():
        # Temporarily update the summary with the bootstrap sample
        curve = fit.calc_pd_curve(phi, df, n_points=20, neighbors=10)
        curve["r"] = np.sqrt(curve.conc_A_rand**2 + curve.conc_B_rand**2)
        ax.plot(curve.r, curve.psi_est, "o", alpha=0.1, color=color)
        a = round(np.cos(phi), 2)
        b = round(np.sin(phi), 2)
        label = f"{fit.drug_A} {a/max(a,b)} : {b/max(a,b)} {fit.drug_B} "
        labels.append(label)
# Create custom legend handles for the three angles
handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor=angle_colors[np.pi / 4],
        markersize=8,
        label=f"$\phi=\pi/4$, {labels[0]}",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor=angle_colors[0],
        markersize=8,
        label=f"$\phi=0$, {labels[1]}",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor=angle_colors[np.pi / 2],
        markersize=8,
        label=f"$\phi=\pi/2$, {labels[2]}",
    ),
]
ax.legend(handles=handles)
ax.set_ylim(-30, 1.5)
plt.show()


import numpy as np
import matplotlib.pyplot as plt

y0 = 10
x1 = 5
phi = [0, np.pi / 16, np.pi / 8, np.pi / 4, np.pi / 2]
m = -y0 / x1
x = y0 / (np.cos(phi) - m)

plt.plot(
    x,
)


x0 = 5
x = np.linspace(0, x0, 100)
y0 = 10
m = -y0 / x0
y = m * x + y0

z, phi = cartesian_to_polar(x, y)

plt.plot(x, y)
plt.scatter(*polar_to_cartesian(z, phi), color="red")


def polar_to_cartesian(r, phi):
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return x, y


def cartesian_to_polar(x, y):
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return r, phi


plt.plot(phi, z)
