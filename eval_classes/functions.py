import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from .bootstrap import BootStrapRamp
from .fit_reps import FitReplicates
from .permutation import PermutationTest
import numpy as np
from icecream import ic
from .solvers import solvers


def call(data, id_value, signal_type, method_name, id_col="drug_id"):
    method = method_name
    method_args = solvers[method_name]
    df = data[(data[id_col] == id_value) & (data.signal_type == signal_type)]
    y = np.array(df["signal"])
    t = np.array(df.treatment_duration)

    if len(np.unique(t)) > 3:
        bootstrap = BootStrapRamp()
        bootstrap.initialize_from_scratch(t, y, method, method_args)
        return bootstrap
    else:
        return False


def collect(data, method, method_args, exclude=[]):
    results = []
    fits = []
    for exp_id in data.exp_id.unique():
        if exp_id not in exclude:
            ic(exp_id)
            for signal_type in ["cfu", "rlu"]:
                ic(signal_type)

                res = call(data, exp_id, signal_type, method, method_args)
                if res:
                    row = res.rate_stats
                    row.update(
                        {"antibiotic": exp_id, "signal_type": signal_type, "obj": res}
                    )
                    results.append(row)
                    fits.append(res)
    return fits, pd.DataFrame().from_records(results)


def t_plot(
    results,
    antibiotic,
    ax=None,
    figsize=(12, 8),
    rlu_color="magenta",
    cfu_color="cyan",
    shape_color="black",
):
    if not ax:
        _, ax = plt.subplots(figsize=figsize)
    sub = results[results.antibiotic == antibiotic]
    rlu_info = sub[sub.signal_type == "rlu"]
    cfu_info = sub[sub.signal_type == "cfu"]
    if len(rlu_info):
        rlu_fit = rlu_info.obj.values[0]
        rlu_fit.plot(ax=ax, color=rlu_color)
    if len(cfu_info):
        cfu_fit = cfu_info.obj.values[0]
        cfu_fit.plot(ax=ax, color=cfu_color)


def get_permutation(results, antibiotic):
    sub = results[results.antibiotic == antibiotic]
    rlu_info = sub[sub.signal_type == "rlu"]
    cfu_info = sub[sub.signal_type == "cfu"]
    if rlu_info.empty or cfu_info.empty:
        print("empty")
        return None
    else:
        rlu_obj = rlu_info.obj.values[0]
        cfu_obj = cfu_info.obj.values[0]
        return PermutationTest(cfu_obj.df.x1.values, rlu_obj.df.x1.values)


def load_data(results, antibiotic):
    sub = results[results.antibiotic == antibiotic]
    rlu_info = sub[sub.signal_type == "rlu"]
    cfu_info = sub[sub.signal_type == "cfu"]

    rlu_fit = rlu_info.obj.values[0] if len(rlu_info) else None
    cfu_fit = cfu_info.obj.values[0] if len(cfu_info) else None

    rlu_x1 = rlu_fit.df["x1"] if rlu_fit else None
    cfu_x1 = cfu_fit.df["x1"] if cfu_fit else None

    return rlu_x1, cfu_x1, rlu_fit, cfu_fit


def plot_stats(ax, fit, position, color="black"):
    if fit:
        ax.scatter(x=[position], y=fit.rate_stats["mean"], color=color, label="Mean")
        ax.errorbar(
            x=[position],
            y=fit.rate_stats["mean"],
            yerr=[
                [fit.rate_stats["mean"] - fit.rate_stats["lower_ci"]],
                [fit.rate_stats["upper_ci"] - fit.rate_stats["mean"]],
            ],
            fmt="none",
            color=color,
            capsize=5,
            label="95% CI",
        )


def violine(ax, combined_df, rlu_fit, cfu_fit, rlu_color="magenta", cfu_color="cyan"):
    sns.violinplot(
        data=combined_df,
        inner=None,
        ax=ax,
        palette={"RLU_x1": rlu_color, "CFU_x1": cfu_color},
    )
    plot_stats(ax, rlu_fit, 0)
    plot_stats(ax, cfu_fit, 1)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["RLU", "CFU"])
    ax.set_xlabel("Signal Type")
    ax.set_ylabel("x1 Value")


def violine_plot(
    results, antibiotic, ax=None, figsize=(12, 8), rlu_color="magenta", cfu_color="cyan"
):
    if not ax:
        _, ax = plt.subplots(figsize=figsize)

    rlu_x1, cfu_x1, rlu_fit, cfu_fit = load_data(results, antibiotic)

    if rlu_x1 is not None and cfu_x1 is not None:
        combined_df = pd.concat(
            [rlu_x1.rename("RLU_x1"), cfu_x1.rename("CFU_x1")], axis=1
        )
        violine(ax, combined_df, rlu_fit, cfu_fit, rlu_color, cfu_color)

    plt.legend()


def get_rep_rates(data, antibiotics, method, method_args):
    results = []
    for antibiotic in antibiotics:
        ic(antibiotic)
        for signal_type in ["cfu", "rlu"]:
            ic(signal_type)
            res = FitReplicates(
                data, antibiotic, signal_type, method, method_args=method_args
            )
            row = {"antibiotic": antibiotic, "signal_type": signal_type, "obj": res}
            results.append(row)
    return pd.DataFrame().from_records(results)


def permute_rates(results, antibiotic, permutation_type="independent"):
    sub = results[results.antibiotic == antibiotic]
    rlu = sub[sub.signal_type == "rlu"]
    rlu_rates = [fits.x[1] for fits in rlu.obj.values[0].fit_results.values()]
    cfu = sub[sub.signal_type == "cfu"]
    cfu_rates = [fits.x[1] for fits in cfu.obj.values[0].fit_results.values()]
    if (len(rlu_rates) > 2) & (len(cfu_rates) > 2):
        permutation = PermutationTest(
            rlu_rates, cfu_rates, permutation_type=permutation_type
        )
        return permutation
    else:
        return None


def plot_rep_rates(
    results,
    antibiotic,
    rlu_color="magenta",
    cfu_color="cyan",
    ax=None,
    figsize=(12, 8),
):
    if not ax:
        _, ax = plt.subplots(figsize=figsize)
    sub = results[results.antibiotic == antibiotic]
    rlu = sub[sub.signal_type == "rlu"]
    cfu = sub[sub.signal_type == "cfu"]
    rlu.obj.values[0].plot(ax=ax, color=rlu_color)
    cfu.obj.values[0].plot(ax=ax, color=cfu_color)


def compare_methods(fits):
    # Assuming 'fits' is your list of fit objects and 'methods' is an array of optimization methods you want to test
    methods = [
        "Nelder-Mead",
        "Powell",
        "CG",
        "BFGS",
        "L-BFGS-B",
        "TNC",
        "COBYLA",
        "SLSQP",
        "trust-constr",
    ]

    # Initialize a list to collect data for each fit and method
    fit_data = []

    # Loop over each fit in your list of fits
    for fit in fits:
        # Loop over each optimization method
        for method in methods:
            # Perform optimization with the current method
            result = fit.optimize(fit.t_data, fit.y_data, method)
            params = result.x
            # Calculate the error using the optimized parameters
            error = fit.error(params, fit.t_data, fit.y_data)

            # Collect information about the fit and optimization result
            fit_info = {
                "Method": method,
                "Params": params,
                "Error": error,
                "Antibiotic": getattr(
                    fit, "antibiotic", "N/A"
                ),  # Assuming 'antibiotic' is an attribute of the fit object
                "FitObject": fit,
            }

            # Append the collected information to our list
            fit_data.append(fit_info)

    # Convert the list to a DataFrame for easy analysis
    df_fits = pd.DataFrame(fit_data)

    # If you need to access mean parameters or any specific attribute, you can now do so from df_fits
    # For example, to get the method with the lowest mean error:
    best_method = df_fits.groupby("Method")["Error"].mean().idxmin()
    print(f"The best fitting method is: {best_method}")
    return df_fits.groupby("Method").mean()


def combine_signals(experiment, img_manager):
    cols = [
        "well",
        "treatment_duration",
        "t",
        "p",
        "antibiotic",
        "signal_type",
        "signal",
    ]
    lum_df = experiment.df.copy()
    lum_df["signal_type"] = "rlu"
    lum_df = lum_df[lum_df.signal.isnull() == False]

    cfu_df = img_manager.df.copy()
    cfu_df["signal_type"] = "cfu"
    cfu_df.rename(columns={"cfu": "signal"}, inplace=True)

    return pd.concat([cfu_df[cols], lum_df[cols]])
