import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


def plot_df2(df, antibiotic, ax1, title=None, lum_signal_type="lum", lum="lum"):
    palette = sns.color_palette("Paired")
    cfu_color = palette[1]
    rlu_color = palette[-5]

    df_a = df[df.antibiotic == antibiotic]

    df_L = df_a[df_a.signal_type == lum_signal_type]
    df_C = df_a[df_a.signal_type == "cfu"]
    ax2 = ax1.twinx()
    _, ax1 = timeplot(
        df_C,
        x_col="t",
        y_col="signal",
        signal_type="cfu",
        hue_col="sec",
        color=cfu_color,
        ax=ax1,
        xlabel="time [h]",
        ylabel="CFU [per ml]",
        lower=10000,
        title=title,
    )
    _, ax2 = timeplot(
        df_L,
        x_col="t",
        y_col=lum,
        signal_type=lum_signal_type,
        hue_col="sec",
        color=rlu_color,
        ax=ax2,
        xlabel="time [h]",
        ylabel="RLU",
        title=title,
    )

    lower_ax1, upper_ax1 = ax1.get_ylim()
    lower_ax2, upper_ax2 = ax2.get_ylim()

    # Apply common limits
    ax1.set_ylim(min(lower_ax1, lower_ax2), max(upper_ax1, upper_ax2))
    ax2.set_ylim(min(lower_ax1, lower_ax2), max(upper_ax1, upper_ax2))


def plot(lum, img_manager, antibiotic, ax1, title=None, lum_label="lum"):
    palette = sns.color_palette("Paired")
    cfu_color = palette[1]
    rlu_color = palette[-5]
    df_L = lum.df[lum.df.antibiotic == antibiotic]
    df_C = img_manager.df[img_manager.df.antibiotic == antibiotic]
    ax2 = ax1.twinx()
    _, ax1 = timeplot(
        df_C,
        x_col="t",
        y_col="signal",
        signal_type="cfu",
        hue_col="sec",
        color=cfu_color,
        ax=ax1,
        xlabel="time [h]",
        ylabel="CFU [per ml]",
        lower=10000,
        title=title,
    )
    _, ax2 = timeplot(
        df_L,
        x_col="t",
        y_col="rlu",
        signal_type=lum_label,
        hue_col="sec",
        color=rlu_color,
        ax=ax2,
        xlabel="time [h]",
        ylabel="RLU",
        title=title,
    )

    lower_ax1, upper_ax1 = ax1.get_ylim()
    lower_ax2, upper_ax2 = ax2.get_ylim()

    # Apply common limits
    ax1.set_ylim(min(lower_ax1, lower_ax2), max(upper_ax1, upper_ax2))
    ax2.set_ylim(min(lower_ax1, lower_ax2), max(upper_ax1, upper_ax2))


def timeplot(
    df,
    x_col,
    y_col,
    signal_type,
    hue_col,
    color,
    ax=None,
    figsize=(12, 8),
    alpha=0.2,
    s=100,
    title=None,
    lower=None,
    upper=None,
    xlabel="t [h]",
    ylabel="Signal",
    legend=True,
):
    df = df[df["signal_type"] == signal_type]

    if not ax:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Check for detection limits
    if lower is not None:
        below_detection_mask = df[y_col] < lower
    else:
        below_detection_mask = pd.Series([False] * len(df), index=df.index)

    if upper is not None:
        above_detection_mask = df[y_col] > upper
    else:
        above_detection_mask = pd.Series([False] * len(df), index=df.index)

    if below_detection_mask.any():
        ax.axhline(lower, color=color, linestyle="--")
        ax.text(
            1.01,
            lower,
            "Lower Detection Limit",
            color=color,
            verticalalignment="bottom",
            horizontalalignment="right",
            transform=ax.get_yaxis_transform(),
        )

    if above_detection_mask.any():
        ax.axhline(upper, color=color, linestyle="--")
        ax.text(
            1.01,
            upper,
            "Upper Detection Limit",
            color=color,
            verticalalignment="top",
            horizontalalignment="right",
            transform=ax.get_yaxis_transform(),
        )

    normal_data = df[~below_detection_mask & ~above_detection_mask]
    below_data = df[below_detection_mask]
    above_data = df[above_detection_mask]

    # Define markers for each well
    markers = [
        "o",
        "s",
        "D",
        "^",
        "v",
        "<",
        ">",
        "p",
        "*",
        "h",
        "H",
        "x",
        "d",
        "|",
        "_",
    ]
    unique_wells = normal_data["well"].unique()
    well_marker_map = {
        well: markers[i % len(markers)] for i, well in enumerate(unique_wells)
    }

    # Plot the data
    for well, marker in well_marker_map.items():
        well_data = normal_data[normal_data["well"] == well]
        ax.scatter(
            well_data[x_col],
            well_data[y_col],
            color=color,
            s=s,
            label=f"{hue_col}: {well}",
            marker=marker,
        )

    ax.scatter(below_data[x_col], below_data[y_col], color=color, s=s, alpha=alpha)
    ax.scatter(above_data[x_col], above_data[y_col], color=color, s=s, alpha=alpha)

    ax.set_yscale("log")
    ax.set_title(title)
    ax.set_ylabel(ylabel, color=color)
    ax.set_xlabel(xlabel)

    # Set y-axis color
    ax.tick_params(axis="y", colors=color)  # Set y-axis ticks color
    # Set y-axis spine color
    ax.yaxis.label.set_color(color)

    # Set y-limits automatically
    lower_ylim = 10 ** np.floor(np.log10(df[y_col][df[y_col] > 0].min()))
    upper_ylim = 10 ** np.ceil(np.log10(df[y_col].max()))
    ax.set_ylim(lower_ylim, upper_ylim)

    return fig, ax
