import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.lines as mlines

paired = sns.color_palette("Paired")
deep = sns.color_palette("deep")
pastel = sns.color_palette("pastel")


def plot_df(
    df,
    antibiotic,
    ax,
    title=None,
    lum_signal_type="lum",
    cfu_color=deep[0],
    rlu_color=paired[-5],
    lum="lum",
    s=100,
    cfu_marker="d",  # New marker for CFU
    rlu_marker="v",  # New marker for RLU
):
    df_a = df[df.antibiotic == antibiotic]

    df_L = df_a[df_a.signal_type == lum_signal_type]
    df_C = df_a[df_a.signal_type == "cfu"]

    _, ax = timeplot(
        df_C,
        x_col="t",
        y_col="signal",
        signal_type="cfu",
        hue_col="sec",
        color=cfu_color,
        marker=cfu_marker,  # Pass CFU marker
        xlabel="time [h]",
        ylabel="signal",
        lower=10000,
        title=title,
        ax=ax,
        s=s,
    )
    _, ax = timeplot(
        df_L,
        x_col="t",
        y_col=lum,
        signal_type=lum_signal_type,
        hue_col="sec",
        color=rlu_color,
        marker=rlu_marker,  # Pass RLU marker
        ax=ax,
        xlabel="time [h]",
        ylabel="signal",
        title=title,
        s=s,
    )

    # Adjust legend marker sizes and include markers
    cfu_legend = mlines.Line2D(
        [],
        [],
        color=cfu_color,
        marker=cfu_marker,  # Use CFU marker
        linestyle="None",
        markersize=np.sqrt(s),
        label="cfu [per ml]",
    )
    rlu_legend = mlines.Line2D(
        [],
        [],
        color=rlu_color,
        marker=rlu_marker,  # Use RLU marker
        linestyle="None",
        markersize=np.sqrt(s),
        label="light intensity [rlu]",
    )

    ax.legend(handles=[cfu_legend, rlu_legend], loc="best", frameon=False)
    ax.set_xlabel("t [h]", labelpad=1)
    ax.set_ylabel("signal", labelpad=1)


def timeplot(
    df,
    x_col,
    y_col,
    signal_type,
    hue_col,
    color,
    ax=None,
    figsize=(12, 8),
    alpha=0.5,
    s=150,
    title=None,
    lower=None,
    upper=None,
    xlabel="t [h]",
    ylabel="Signal",
    marker="o",  # Default marker
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
        ax.axhline(lower, color=color, linestyle=":", alpha=alpha)
        ax.text(
            0.99,
            lower,
            "cfu detection limit",
            alpha=alpha,
            color=color,
            verticalalignment="bottom",
            horizontalalignment="right",
            transform=ax.get_yaxis_transform(),
        )

    if above_detection_mask.any():
        ax.axhline(upper, color=color, linestyle=":")
        ax.text(
            0.99,
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
        ".",
        ".",
        ".",
        ".",
        ".",
        ".",
        ".",
    ]
    unique_wells = normal_data["well"].unique()
    well_marker_map = {
        well: markers[i % len(markers)] for i, well in enumerate(unique_wells)
    }

    ax.scatter(
        normal_data[x_col],
        normal_data[y_col],
        color=color,
        s=s,
        marker=marker,  # Use passed marker
    )

    ax.scatter(
        below_data[x_col],
        below_data[y_col],
        color=color,
        s=s,
        alpha=alpha,
        marker=marker,
    )
    ax.scatter(
        above_data[x_col],
        above_data[y_col],
        color=color,
        s=s,
        alpha=alpha,
        marker=marker,
    )

    ax.set_yscale("log")
    ax.set_title(title)
    ax.set_ylabel(ylabel, labelpad=1)  # Shift ylabel closer
    ax.set_xlabel(xlabel, labelpad=1)  # Shift xlabel closer

    # Set y-limits automatically
    ax.set_ylim(10**2, 10**9)

    return fig, ax
