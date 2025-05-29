import matplotlib.pyplot as plt
import os
from general_classes import TimeLog, PathManager
import os
import sys
import numpy as np
from figures import Figure, styles

figure = Figure(styles["presentation"])


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
summaries = {}
combination_objs = {}

palette = {"active": (1, 1, 0.7), "ended": (0.6, 0.6, 0.6), "fit": (1, 0.7, 0.7)}

for i in [5, 7, 8, 13]:  # 5, 7, 13
    key = keys[i]
    comb_dict = combinations[key]
    comb_obj = FitRates(comb_dict, path_manager)
    summary = comb_obj.load_summary()
    summaries[key] = summary
    combination_objs[key] = comb_obj


drug = "chloramphenicol_tetracycline"
obj = combination_objs[drug]

idx_pair = (0, 0)
fig, ax = figure.create_figure_with_style((12, 8))
obj.plot_concentration_idx_pair(idx_pair, ax=ax, palette=palette)
obj.plot_fit_idx_pair(idx_pair, ax=ax, palette=palette)
plt.show()


def mk_idx_pair_figure(drug, idx_pair):
    obj = combination_objs[drug]
    _, ax = figure.create_figure_with_style((12, 8))
    obj.plot_concentration_idx_pair(idx_pair, ax=ax, palette=palette, fontsize=20)
    obj.plot_fit_idx_pair(idx_pair, ax=ax, palette=palette, fontsize=20)
    figure.save_figure(f"{drug}_{idx_pair}", folderpath="figures")
    plt.show()


mk_idx_pair_figure("chloramphenicol_tetracycline", (0, 0))
mk_idx_pair_figure("fosfomycin_tetracycline", (11, 9))
mk_idx_pair_figure("chloramphenicol_colistin", (0, 11))
mk_idx_pair_figure("chloramphenicol_polymyxinB", (0, 11))


def mk_figures(drug, vmin, vmax):
    obj = combination_objs[drug]

    _, ax = figure.create_figure_with_style((24, 8))
    obj.median_heatmap(save_local=True, ax=ax)
    figure.save_figure(f"{drug}_median", folderpath="figures")
    plt.show()

    _, ax = figure.create_figure_with_style((24, 8))
    obj.plot_plate("I", cbar_shrink=0.7, vmin=vmin, vmax=vmax, ax=ax, fontsize=30)
    figure.save_figure(f"{drug}_I", folderpath="figures")
    plt.show()

    _, ax = figure.create_figure_with_style((24, 8))
    obj.plot_plate("II", cbar_shrink=0.7, vmin=vmin, vmax=vmax, ax=ax, fontsize=30)
    figure.save_figure(f"{drug}_II", folderpath="figures")
    plt.show()

    _, ax = figure.create_figure_with_style((12, 12))
    obj.plot_median_contour(fontsize=30, ax=ax, sigma=100)
    figure.save_figure(f"{drug}_median_contour", folderpath="figures")
    plt.show()


mk_figures("chloramphenicol_tetracycline", vmin=-1, vmax=2)
mk_figures("chloramphenicol_colistin", vmin=-20, vmax=2)
mk_figures("chloramphenicol_polymyxinB", vmin=-50, vmax=2)
mk_figures("fosfomycin_tetracycline", vmin=-5, vmax=2)

obj = combination_objs["chloramphenicol_polymyxinB"]
idx_pair = (0, 0)
fig, ax = obj.plot_concentration_idx_pair(idx_pair)
obj.plot_fit_idx_pair(idx_pair, ax=ax)
plt.show()


idx_pair = (0, 9)
obj = combination_objs["chloramphenicol_colistin"]
fig, ax = obj.plot_concentration_idx_pair(idx_pair)
obj.plot_fit_idx_pair(idx_pair, ax=ax)
plt.show()


idx_pair = (11, 9)
obj = combination_objs["fosfomycin_tetracycline"]
fig, ax = obj.plot_concentration_idx_pair(idx_pair)
obj.plot_fit_idx_pair(idx_pair, ax=ax)
plt.show()


#######################################################
######### PD CURVES
phi = np.pi / 4
palette = {
    "measured": (0.7, 0.7, 0.7),  # light grey
    "fit": (0.975, 0.965, 0.94),  # pearl white
    "loewe": (0, 1, 1.0),  # cyan (soft / light blue)
    "bliss": (1.0, 0.4, 0.38),  # pastel red
}


combi = "fosfomycin_tetracycline"
fit = PdFitter(summaries[combi], palette=palette)
_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(phi=phi, loewe=True, bliss=True, ax=ax, linewidth=4, s=200)
figure.save_figure(combi, folderpath="figures")
plt.show()


_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=False,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi + "_loewe", folderpath="figures")
plt.show()


combi = "chloramphenicol_colistin"
fit = PdFitter(summaries[combi], palette=palette)
_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=True,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi, folderpath="figures")
plt.show()


_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=False,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi + "_loewe", folderpath="figures")
plt.show()


combi = "chloramphenicol_polymyxinB"
fit = PdFitter(summaries[combi], palette=palette)
_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=True,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi, folderpath="figures")
plt.show()

_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=False,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi + "_loewe", folderpath="figures")
plt.show()


combi = "chloramphenicol_tetracycline"
fit = PdFitter(summaries[combi], palette=palette)
_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=True,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi, folderpath="figures")
plt.show()

_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=False,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi + "_loewe", folderpath="figures")
plt.show()


_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=False,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi + "_loewe", folderpath="figures")
plt.show()

_, ax = figure.create_figure_with_style()
fit.plot_mono("A", ax=ax)
figure.save_figure(fit.drug_A + "_pdA", folderpath="figures")
plt.show()

_, ax = figure.create_figure_with_style()
fit.plot_mono("B", ax=ax)
figure.save_figure(fit.drug_B + "_pdB", folderpath="figures")
plt.show()


_, ax = figure.create_figure_with_style()
fit.plot_combi_pd(
    phi=phi,
    loewe=False,
    bliss=False,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi + "_data", folderpath="figures")
plt.show()

fit.plot_combi_pd(
    phi=phi,
    loewe=False,
    bliss=True,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi + "_bliss", folderpath="figures")

fit.plot_combi_pd(
    phi=phi,
    loewe=True,
    bliss=False,
    ax=ax,
    linewidth=4,
    s=200,
)
figure.save_figure(combi + "_loewe", folderpath="figures")
