from scipy.stats import permutation_test
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import permutation_test
import numpy as np
import matplotlib.pyplot as plt


class PermutationTest:
    def __init__(
        self,
        X1,
        X2,
        n_resamples=9999,
        alternative="two-sided",
        stat_type="mean",
        permutation_type="independent",
    ):
        self.X1 = X1
        self.X2 = X2
        if stat_type == "mean":
            statistic = lambda x, y: np.mean(x) - np.mean(y)
        elif stat_type == "median":
            statistic = lambda x, y: np.median(x) - np.median(y)
        else:
            raise ValueError("stat_type must be 'mean' or 'median'")

        self.result = permutation_test(
            (X1, X2),
            statistic=statistic,
            permutation_type=permutation_type,
            n_resamples=n_resamples,
            alternative=alternative,
        )
        self.pvalue = self.result.pvalue
        self.null_distribution = self.result.null_distribution

    def plot(self, ax=None, figsize=(12, 8), bins=10):
        if not ax:
            _, ax = plt.subplots(figsize=figsize)
        significance_level = 0.05
        color = "green" if self.pvalue > significance_level else "red"
        label = f"Null distribution"

        ax.hist(self.null_distribution, bins=bins, alpha=0.7, label=label)
        ax.axvline(
            self.result.statistic,
            color="black",
            linestyle="dashed",
            linewidth=2,
            label="Observed difference",
        )
        # ax.legend()
        ax.set_xlabel("Difference in means")
        ax.set_ylabel("Frequency")
        ax.set_title(f"Permutation - p-value: {self.pvalue:.4f}", color=color)
