import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from icecream import ic
import pickle


class BootStrapRamp:
    def __init__(self):
        """Initialize an empty BootStrapRamp instance."""
        pass

    def initialize_from_scratch(
        self,
        t_data,
        y_data,
        method,
        cfu_color="green",
        rlu_color="red",
        method_args=None,
        dt_min=2,
    ):
        """Initialize the BootStrapRamp instance with data and computations."""
        self.t_data = t_data
        self.y_data = y_data
        self.cfu_color = cfu_color
        self.rlu_color = rlu_color
        self.method = method
        self.method_args = method_args if method_args is not None else {}
        self.df = self.sample_data()
        self.initial_guess = self.calculate_initial_guess(t_data, y_data, dt_min)
        self.bounds = self.calculate_bounds(t_data, y_data, dt_min)
        self.optimized_bootstrapped(self.df)
        self.rate_stats = self.get_stats()
        self.t_fit, self.y_fit, self.mean_params = self.mean_fit()

    def error(self, params, t_data, y_data):
        y_sim = self.ramp_function(t_data, params)
        e1 = np.sum((np.log(y_data + 1) - np.log(y_sim + 1)) ** 2)
        mask = (params[2] <= t_data) & (t_data <= params[3])
        e2 = np.sum(((np.log(y_data[mask]) - np.log(y_sim[mask])) / sum(mask)) ** 2)
        return e1

    def ramp_function(self, T, params):
        y_sim = [self.ramp_function_t(t, params) for t in T]
        return np.array(y_sim)

    @staticmethod
    def ramp_function_t(t, params):
        y_0, rate, t0, dt = params
        y_p = y_0 * np.exp(-rate * t0)
        te = t0 + dt
        if t < t0:
            y = y_p * np.exp(rate * t0)
        elif t > te:
            y = y_p * np.exp(rate * te)
        else:
            y = y_p * np.exp(rate * t)
        return y

    def calculate_initial_guess(self, t_data, y_data, dt_min):
        t0 = np.min(t_data)
        t_end = np.max(t_data)
        Y0 = y_data[t_data == t0]
        Y_end = y_data[t_data == t_end]

        # Initial y0 guess
        y0_initial_guess = np.exp(np.mean(np.log(Y0)))

        # Improve initial rate calculation
        mean_start = np.mean(Y0)
        mean_end = np.mean(Y_end)
        rate_initial_guess = (np.log(mean_end) - np.log(mean_start)) / (t_end - t0)

        dt = dt_min * 1.25
        initial_guess = [y0_initial_guess, rate_initial_guess, t_end / 3, dt]

        return initial_guess

    def calculate_bounds(self, t_data, y_data, dt_min):
        Y0 = y_data[t_data == np.min(t_data)]
        bounds = [
            (min(Y0) / 4, max(Y0) * 4),
            (-20, 5),  # Assuming rate is negative; adjust based on your data
            (min(t_data), max(t_data)),
            (
                dt_min,
                np.max(t_data) - np.min(t_data),
            ),  # dt bound based on the total duration
        ]
        return bounds

    def optimize(self, t_data, y_data, method):

        result = minimize(
            self.error,
            self.initial_guess,
            args=(t_data, y_data),
            method=method,
            options={**self.method_args},
            bounds=self.bounds,
        )
        return result

    def sample_data(self, num_samples=100):
        bootstrapped_data = []
        for i in range(num_samples):
            li = 0
            while li < 3:
                indices = np.random.choice(
                    len(self.t_data), size=len(self.t_data), replace=True
                )
                li = len(np.unique(indices))

            name = "set" + str(i)
            t_bootstrap = self.t_data[indices]
            y_bootstrap = self.y_data[indices]
            bootstrapped_data.append({"set": name, "t": t_bootstrap, "y": y_bootstrap})
        return pd.DataFrame().from_records(bootstrapped_data).set_index("set")

    def optimized_bootstrapped(self, sampled_df):
        for i, row in sampled_df.iterrows():
            params = self.optimize(row.t, row.y, self.method)
            x0 = params.x[0]
            x1 = params.x[1]
            x2 = params.x[2]
            x3 = params.x[3]
            sampled_df.loc[i, ["x0", "x1", "x2", "x3", "params"]] = [
                x0,
                x1,
                x2,
                x3,
                params,
            ]

    def get_stats(self, confidence_level=95):
        rate_stats = {
            "mean": self.df.x1.mean(),
            "var": self.df.x1.var(),
            "lower_ci": np.percentile(self.df.x1, (100 - confidence_level) / 2, axis=0),
            "upper_ci": np.percentile(
                self.df.x1, confidence_level + (100 - confidence_level) / 2, axis=0
            ),
        }
        return rate_stats

    def mean_fit(self, n=100):
        mean_params = [
            self.df.x0.mean(),
            self.df.x1.mean(),
            self.df.x2.mean(),
            self.df.x3.mean(),
        ]
        t_fit = np.linspace(min(self.t_data), max(self.t_data), n)
        y_fit = self.ramp_function(t_fit, mean_params)
        return t_fit, y_fit, mean_params

    def plot(self, ax=None, figsize=(12, 8), color="magenta"):
        if not ax:
            _, ax = plt.subplots(figsize=figsize)
        ax.scatter(self.t_data, self.y_data, color=color)
        ax.plot(self.t_fit, self.y_fit, color=color)
        ax.set_xlabel("t")
        ax.set_ylabel("y")
        ax.set_yscale("log")

    def save_to_file(self, filename):
        """Save the instance to a file using pickle."""
        with open(filename, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load_from_file(filename):
        """Load an instance from a file using pickle."""
        with open(filename, "rb") as f:
            return pickle.load(f)
