import matplotlib.pyplot as plt
import seaborn as sns

pastel = sns.color_palette("pastel")

# Define color schemes for "presentation" and "paper" modes
color_schemes = {
    "paper": {
        "rlu_color": "#FFA500",  # bright orange
        "cfu_color": "#000080",  # dark navy blue
    },
    "presentation": {
        "rlu_color": pastel[-2],  # pastel yellow
        "cfu_color": pastel[-1],  # pastel blue
    },
}


def cfu_rlu_timelplot(
    df,
    antibiotic,
    ax=None,
    figsize=(12, 8),
    title=None,
    color_mode="paper",
    legend=True,
):
    if color_mode not in color_schemes:
        raise ValueError(f"color_mode must be one of {list(color_schemes.keys())}")

    colors = color_schemes[color_mode]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    sub = df[df.antibiotic == antibiotic]
    ax.set_yscale("log")

    rlu = sub[sub.signal_type == "rlu"]
    ax.scatter(
        rlu.treatment_duration, rlu.signal, color=colors["rlu_color"], label="RLU"
    )

    cfu = sub[sub.signal_type == "cfu"]
    ax.scatter(
        cfu.treatment_duration, cfu.signal, color=colors["cfu_color"], label="CFU"
    )

    ax.set_ylabel("Signal")
    ax.set_xlabel("Time [h]")
    ax.set_title(title)

    if legend:
        ax.legend()

    return ax
