import pandas as pd
import numpy as np
import json
from .morph_correction import morph_correction

morph_path = "/Users/malte/ETH-Documents/Luminescence_Paper/morphology/summary_df.csv"
morph_summary_r = (
    pd.read_csv(morph_path)[["antibiotic", "replicate", "volume"]]
    .groupby(["antibiotic", "replicate"])
    .mean()
)
morph_summary = morph_summary_r.groupby("antibiotic").mean()


class CompareSlopes:
    def __init__(
        self,
        info,
        pm,
        significant_filamentation=False,
        bootstrap_path="/Users/malte/ETH-Documents/Luminescence_Paper/comparison_code/bootstrapped/",
    ):
        name_add = info["name_add"]
        self.bootstrap_path = bootstrap_path
        self.pm = pm
        self.info = info
        self.antibiotic = info["antibiotic"]
        self.df, self.df_below = self._load_filtered_data()
        self.t_cfu, self.y_cfu, self.t_lum, self.y_lum = self._extract_signals(
            info["lum_label"]
        )

        self.cfu_fits = self._load_fit_results("_bootstraped_cfu_slopes", name_add)
        self.lum_fits = self._load_fit_results("_bootstraped_lum_slopes", name_add)

        lum_summary, self.lum_slopes = self.summarize(self.lum_fits)
        cfu_summary, self.cfu_slopes = self.summarize(self.cfu_fits)

        self.slope_df = self._mk_df()
        self.df = self.slope_df.copy()
        self.df["antibiotic"] = info["antibiotic"]
        self.summary = self._compute_summary(lum_summary, cfu_summary)

        if bool(info.get("morph_name", None)) & significant_filamentation:
            self._process_morph_correction(name_add)
        else:
            self.morph = False

    def _load_filtered_data(self):
        file_path = self.pm.file_path(
            self.info["label"] + self.info["name_add"] + "_df_filtered.csv",
            folder="eval_code",
        )
        df = pd.read_csv(file_path)
        return (
            df[(df.antibiotic == self.antibiotic) & df.include],
            df[(df.antibiotic == self.antibiotic) & df.below_limit],
        )

    def _extract_signals(self, lum_label):
        df = self.df
        t_cfu = df[df.signal_type == "cfu"].t.to_numpy()
        y_cfu = df[df.signal_type == "cfu"].signal.to_numpy()
        t_lum = df[df.signal_type == lum_label].t.to_numpy()
        y_lum = df[df.signal_type == lum_label].signal.to_numpy()
        return t_cfu, y_cfu, t_lum, y_lum

    def _load_fit_results(self, suffix, name_add):
        file_name = f"{self.antibiotic}{suffix}{name_add}.json"
        return self.load_json(self.bootstrap_path + file_name)

    def _compute_summary(self, lum_summary, cfu_summary):
        return {
            "lum": lum_summary,
            "cfu": cfu_summary,
            "significance": self.significance(self.lum_slopes),
            "abs_diff": abs(lum_summary["mean_slope"] - cfu_summary["mean_slope"]),
        }

    def _process_morph_correction(self, name_add):
        self.morph = True
        morph_file = (
            f"{self.antibiotic}_bootstraped_lum_slopes_morph_corrected{name_add}.json"
        )
        self.lum_fits_corrected = self.load_json(self.bootstrap_path + morph_file)
        self.y_lum_corrected = morph_correction(
            self.t_lum, self.y_lum, self.info["morph_name"], morph_summary
        )
        lum_corr_summary, self.lum_corr_slopes = self.summarize(self.lum_fits_corrected)
        self.summary.update(
            {
                "lum_corr": lum_corr_summary,
                "significance_corr": self.significance(self.lum_corr_slopes),
                "abs_diff": abs(
                    lum_corr_summary["mean_slope"] - self.summary["cfu"]["mean_slope"]
                ),
            }
        )
        corrected_data = pd.DataFrame(
            {
                "slope": self.lum_corr_slopes,
                "signal_type": ["lum_corr"] * len(self.lum_corr_slopes),
            }
        )
        self.slope_df = pd.concat([self.slope_df, corrected_data])

    def _mk_df(self):
        data = {
            "slope": np.concatenate([self.lum_slopes, self.cfu_slopes]),
            "signal_type": (
                ["lum"] * len(self.lum_slopes) + ["cfu"] * len(self.cfu_slopes)
            ),
        }
        return pd.DataFrame(data)

    def load_json(self, file_name):
        """
        Load a JSON file and return its content as a dictionary.

        Parameters:
        file_name (str): Path to the JSON file.

        Returns:
        dict: Dictionary loaded from the JSON file.
        """
        try:
            with open(file_name, "r") as f:
                data = json.load(f)
            # Convert lists back to NumPy arrays where applicable
            for key, value in data.items():
                for k, v in value.items():
                    if isinstance(v, list):
                        data[key][k] = np.array(v)
            return data
        except FileNotFoundError:
            raise FileNotFoundError(f"File {file_name} not found.")
        except json.JSONDecodeError:
            raise ValueError(f"File {file_name} is not a valid JSON file.")

    def summarize(self, fits):
        slopes = []
        intercepts = []
        for dataset in fits.values():
            slopes.append(dataset["slope"])
            intercepts.append(dataset["intercept"])
        return (
            {
                "slopes": slopes,
                "intercepts": intercepts,
                "mean_slope": np.mean(slopes),
                "mean_intercept": np.mean(intercepts),
                "pi": self.compute_pi(slopes),
            },
            slopes,
        )

    def significance(self, lum_slopes):
        lum_mean = np.mean(lum_slopes)
        cfu_mean = np.mean(self.cfu_slopes)

        lum_pi = self.compute_pi(lum_slopes)
        cfu_pi = self.compute_pi(self.cfu_slopes)

        if (cfu_pi[0] <= lum_mean <= cfu_pi[1]) and (
            lum_pi[0] <= cfu_mean <= lum_pi[1]
        ):
            return "n.s."  # Mean of one is within the PI of the other
        else:
            return "*"  # PIs overlap but means don't

    @staticmethod
    def compute_pi(data, confidence=0.95):
        """
        Compute the prediction interval for the given data.

        Parameters:
        data (array-like): The data for which to compute the interval.
        confidence (float): Confidence level for the interval (default: 0.95).

        Returns:
        tuple: Lower and upper bounds of the prediction interval.
        """
        lower = np.percentile(data, (1 - confidence) / 2 * 100)
        upper = np.percentile(data, (1 + confidence) / 2 * 100)
        return lower, upper

    # Additional methods (load_lum_tail, load_json, summarize, significance, compute_pi, etc.) remain unchanged.
