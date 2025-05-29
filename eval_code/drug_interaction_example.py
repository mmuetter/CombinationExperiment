import numpy as np
import matplotlib.pyplot as plt
from figures import Figure, styles
from matplotlib.lines import Line2D

figure = Figure(styles["presentation"])


class PDIsoboles:
    def __init__(self, palette=None):
        default_palette = {"loewe": "C1", "bliss": "C2"}
        self.palette = palette or default_palette
        self.params = {}

    def define_pharmacodynamic_parameters(self, drug, zmic, kappa, psi_min, psi_max):
        self.params[drug] = dict(
            zmic=zmic, kappa=kappa, psi_min=psi_min, psi_max=psi_max
        )

    @staticmethod
    def tau(d, zmic, kappa, psi_min, psi_max):
        num = (psi_max - psi_min) * (d / zmic) ** kappa
        den = (d / zmic) ** kappa - (psi_min / psi_max)
        return num / den

    @staticmethod
    def inv_dose(tau, zmic, kappa, psi_min, psi_max):
        frac = (tau * psi_min / psi_max) / (tau - (psi_max - psi_min))
        return zmic * frac ** (1 / kappa)

    def calc_loewe_isobole(self, tau, num=200):
        pA = self.params["A"]
        pB = self.params["B"]
        DA = self.inv_dose(tau, **pA)
        DB = self.inv_dose(tau, **pB)
        cA = np.linspace(0, DA, num)
        cB = DB * (1 - cA / DA)
        return cA, cB

    def calc_bliss_isobole(self, tau, num=200):
        pA = self.params["A"]
        pB = self.params["B"]
        DA = self.inv_dose(tau, **pA)
        cA = np.linspace(0, DA, num)
        tauA = self.tau(cA, **pA)
        tauB = tau - tauA
        cB = np.array([self.inv_dose(tb, **pB) if tb > 0 else 0 for tb in tauB])
        return cA, cB

    def plot_isoboles(
        self, taus, method="both", logscale=False, fontsize=None, ax=None
    ):
        """
        Plot isoboles for given list of tau values.
        method: 'loewe', 'bliss', or 'both'
        logscale: use log-log axes
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        methods = {"loewe": self.calc_loewe_isobole, "bliss": self.calc_bliss_isobole}
        keys = ["loewe", "bliss"] if method == "both" else [method]

        for key in keys:
            calc = methods[key]
            color = self.palette[key]
            for idx, tau in enumerate(taus):
                try:
                    cA, cB = calc(tau)
                except Exception:
                    continue
                ax.plot(cA, cB, color=color, alpha=1 / (1 + 0.5 * idx))

        ax.set_xlabel(r"$c_A$", fontsize=fontsize)
        ax.set_ylabel(r"$c_B$", fontsize=fontsize)
        if logscale:
            ax.set_xscale("log")
            ax.set_yscale("log")

        # custom legend for models only, alpha=1
        handles = [
            Line2D([0], [0], color=self.palette["loewe"], lw=2, alpha=1, label="Loewe"),
            Line2D([0], [0], color=self.palette["bliss"], lw=2, alpha=1, label="Bliss"),
        ]
        ax.legend(handles=handles, fontsize=fontsize)

        if fontsize:
            ax.tick_params(labelsize=fontsize)
        plt.tight_layout()
        return ax


palette = {"loewe": "cyan", "bliss": (1, 0.6, 0.6)}
# Example usage of PDIsoboles class to plot isoboles for two drugs A and B
iso = PDIsoboles(palette=palette)
iso.define_pharmacodynamic_parameters("A", zmic=1.0, kappa=1.5, psi_min=-10, psi_max=2)
iso.define_pharmacodynamic_parameters("B", zmic=1, kappa=1.5, psi_min=-10, psi_max=2)


iso.plot_isoboles(taus=[8, 4, 2, 1], method="both", logscale=True, fontsize=14)
plt.show()
figure.save_figure("isoboles_log")


iso.plot_isoboles(taus=[8, 4, 2, 1], method="both", logscale=False, fontsize=14)
plt.show()
figure.save_figure("isoboles_linear")
