import numpy as np
import pandas as pd
from .morph_correction import morph_correction
import json

morph_path = "/Users/malte/ETH-Documents/Luminescence_Paper/morphology/summary_df.csv"
morph_summary_r = (
    pd.read_csv(morph_path)[["antibiotic", "replicate", "volume"]]
    .groupby(["antibiotic", "replicate"])
    .mean()
)
morph_summary = morph_summary_r.groupby("antibiotic").mean()


class BootstrapSlope:
    def __init__(
        self,
        info,
        n,
        pm,
        fit_method,
        name_add="",
        cfu_label="cfu",
        windowsize=None,
        m_min=3,
        save_folder="/Users/malte/ETH-Documents/Luminescence_Paper/comparison_code/bootstrapped/",
        significant_filamentation: bool = False,
    ):
        self.windowsize = windowsize
        if windowsize:
            self.m_min = windowsize
        else:
            self.m_min = m_min
        self.fit_method = fit_method

        df = pd.read_csv(
            pm.file_path(
                info["label"] + info["name_add"] + "_df_filtered.csv",
                folder="eval_code",
            )
        )
        self.df = df[df.include]

        self.pm = pm
        self.info = info
        self.lum_label = info["lum_label"]

        self.t_cfu = self.df[self.df.signal_type == cfu_label].t.to_numpy()
        self.t_lum = self.df[self.df.signal_type == info["lum_label"]].t.to_numpy()
        self.y_cfu = self.df[self.df.signal_type == cfu_label].signal.to_numpy()
        self.y_lum = self.df[self.df.signal_type == info["lum_label"]].signal.to_numpy()

        self.i_cfu = self.bootstrap_indices(self.t_cfu, n)
        self.i_lum = self.bootstrap_indices(self.t_lum, n)
        self.cfu_fits = self.fit(self.i_cfu, self.t_cfu, self.y_cfu)
        self.lum_fits = self.fit(self.i_lum, self.t_lum, self.y_lum)

        antibiotic = info["antibiotic"]
        self.save_bootstrap_results(
            self.cfu_fits,
            antibiotic + "_bootstraped_cfu_slopes" + name_add,
            save_folder,
        )
        self.save_bootstrap_results(
            self.lum_fits,
            antibiotic + "_bootstraped_lum_slopes" + name_add,
            save_folder,
        )

        if bool(info.get("morph_name", None)) & significant_filamentation:
            self.y_corr_fits = morph_correction(
                self.t_lum, self.y_lum, info.get("morph_name", None), morph_summary
            )
            self.lum_fits_corr = self.fit(self.i_lum, self.t_lum, self.y_corr_fits)
            self.save_bootstrap_results(
                self.lum_fits_corr,
                antibiotic + "_bootstraped_lum_slopes_morph_corrected" + name_add,
                save_folder,
            )

    def bootstrap_indices(self, t, n):
        """
        Generate `n` bootstrap datasets with at least two unique timepoints.
        """
        indices = []
        for _ in range(n):
            while True:
                sampled_indices = np.random.choice(len(t), size=len(t), replace=True)
                if (
                    len(np.unique(t[sampled_indices])) >= self.m_min
                ):  # Ensure at least 2 unique timepoints
                    indices.append(sampled_indices)
                    break
        return indices

    def fit(self, bootstraped_sets, t, y):
        fits = {}

        for index, indices in enumerate(bootstraped_sets):
            t_sampled = t[indices]
            y_sampled = y[indices]
            if self.windowsize:
                fit = self.fit_method(t_sampled, y_sampled, windowsize=self.windowsize)
            else:
                fit = self.fit_method(t_sampled, y_sampled)
            fits.update({index: fit.summary})

        return fits

    def save_bootstrap_results(self, fits, name, save_folder):
        serializable_fits = {
            key: {
                k: v.tolist() if isinstance(v, np.ndarray) else v
                for k, v in value.items()
            }
            for key, value in fits.items()
        }

        # Save the serialized dictionary to a JSON file
        with open(f"{save_folder}{name}.json", "w") as f:
            json.dump(serializable_fits, f, indent=4)
