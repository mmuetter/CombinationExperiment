import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from icecream import ic


class FitReplicates:
    def __init__(
        self,
        df,
        antibiotic,
        y_label,
        method,
        t_label="treatment_duration",
        rep_label="rep",
        cfu_color="green",
        rlu_color="red",
        method_args=None,
        dt_min=2,
    ):
        self.df = df[df.antibiotic == antibiotic]
        self.cfu_color = cfu_color
        self.rlu_color = rlu_color
        self.method = method

        self.t_label, self.y_label, self.rep_label, self.antibiotic, self.dt_min = (
            t_label,
            y_label,
            rep_label,
            antibiotic,
            dt_min,
        )
        self.method_args = method_args if method_args is not None else {}
        self.initial_guess = self.calculate_initial_guess()
        self.bounds = self.calculate_bounds()
        self.calculate_bounds()
        self.fit_results = self.fit_rates()

    def calculate_initial_guess(self):
        t_data = self.df[self.t_label]
        t0 = np.min(t_data)
        t_end = np.max(t_data)

        y_data = self.df[self.y_label]
        Y0 = y_data[t_data == t0]
        Y_end = y_data[t_data == t_end]

        # Initial y0 guess
        y0_initial_guess = np.exp(np.mean(np.log(Y0)))

        # Improve initial rate calculation
        mean_start = np.mean(Y0)
        mean_end = np.mean(Y_end)
        rate_initial_guess = (np.log(mean_end) - np.log(mean_start)) / (t_end - t0)

        dt = self.dt_min * 1.25
        initial_guess = [y0_initial_guess, rate_initial_guess, t_end / 3, dt]

        return initial_guess

    def calculate_bounds(self):
        t_data = self.df[self.t_label]
        y_data = self.df[self.y_label]

        Y0 = y_data[t_data == np.min(t_data)]
        bounds = [
            (min(Y0) / 4, max(Y0) * 4),
            (-20, 5),  # Assuming rate is negative; adjust based on your data
            (min(t_data), max(t_data)),
            (
                self.dt_min,
                np.max(t_data) - np.min(t_data),
            ),  # dt bound based on the total duration
        ]
        return bounds

    def fit_rates(self):
        rates = {}
        for rep in self.df[self.rep_label].unique():
            sub = self.df[self.df[self.rep_label] == rep]
            if len(sub) > 3:
                result = self.optimize(sub[self.t_label], sub[self.y_label])
                rates.update({rep: result})
        return rates

    def optimize(self, t_data, y_data):
        result = minimize(
            self.error,
            self.initial_guess,
            args=(t_data, y_data),
            method=self.method,
            options={**self.method_args},
            bounds=self.bounds,
        )
        return result

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

    def plot(self, ax=None, figsize=(12, 8), color="magenta", n=100):
        if not ax:
            _, ax = plt.subplots(figsize=figsize)

        ax.scatter(self.df[self.t_label], self.df[self.y_label], color=color)
        self.t = np.linspace(min(self.df[self.t_label]), max(self.df[self.t_label]), n)
        for rep, result in self.fit_results.items():
            y = self.ramp_function(self.t, result.x)
            ax.plot(self.t, y, label="rep" + str(rep), color=color)
        ax.set_xlabel("t")
        ax.set_ylabel("y")
        ax.set_yscale("log")
