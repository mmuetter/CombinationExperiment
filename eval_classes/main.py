from eval_val_curves import ImManager, Experiment
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def indexing(x):
    X = x.assay_id.split("_")
    X.append(x.well)
    return "_".join(X)


def assign_control_wells(well):
    return not well in p0_antibiotics.keys()


experiment = Experiment()
experiment.df["control"] = experiment.df.antibiotic.isnull()
experiment.control = experiment.df[experiment.df.control]
experiment.df = experiment.df[experiment.df.control == False]
experiment.df["index"] = experiment.df["file_name"] + "_" + experiment.df.well
experiment.df.set_index("index", inplace=True)

img_manager = ImManager(experiment, "agar_t5_p0_d1_n1.png")
img_manager.cut_all_sections()
img_manager.select_sections()
img_manager.evaluate_colonies()
img_manager.get_cfu()
img_manager.df["index"] = img_manager.df.apply(lambda x: indexing(x), axis=1)
img_manager.df.set_index("index", inplace=True)


cor = img_manager.reeval_section("t1", "p1", "sec5")
img_manager.get_cfu()

merged_df = pd.merge(
    experiment.df,
    img_manager.df,
)
merged_df.to_csv("merged_df.csv")

merged_df = merged_df.rename(columns={"signal": "rlu"})
for antibiotic in merged_df.antibiotic.unique():
    ab_df = merged_df[merged_df.antibiotic == antibiotic]
    for well in ab_df.well.unique():
        sub_df = ab_df[(ab_df.well == well)]
        sub_t0 = sub_df[sub_df.t == "t0"]
        merged_df.loc[sub_df.index, "rlu_normed"] = sub_df["rlu"] / sub_t0.rlu.values[0]
        merged_df.loc[sub_df.index, "cfu_normed"] = sub_df["cfu"] / sub_t0.cfu.values[0]

id_vars = [col for col in merged_df.columns if col not in ["cfu_normed", "rlu_normed"]]
long_df = pd.melt(
    merged_df,
    id_vars=id_vars,
    value_vars=["cfu_normed", "rlu_normed"],
    var_name="signal_type",
    value_name="signal",
)

g = sns.relplot(
    data=long_df,
    x="treatment_duration",
    y="signal",
    hue="signal_type",
    col="antibiotic",
    col_wrap=2,
    style="sec",
)
g.set(yscale="log")
plt.savefig(img_manager.file_path("normed_signal.png", "analysis"))
plt.show()

# t, p, sec = "t0", "p0", "sec0"
