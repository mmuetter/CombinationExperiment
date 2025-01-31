import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class WindowFit:
    def __init__(self, t, y, windowsize=4):
        self.t = t
        self.y = y
        self.windowsize = windowsize
        self.popt, self.t_fit, self.y_fit = self.find_min_slope()
        self.summary = self.summarize()

    def exp_decline(self, t, a, b):
        return a + b * t

    def initial_guess(self, t, y):
        t_min, t_max = np.min(t), np.max(t)
        y_min = np.mean(y[t == t_min])
        y_max = np.mean(y[t == t_max])
        log_y_min, log_y_max = np.log(y_min), np.log(y_max)
        b_initial = (log_y_max - log_y_min) / (t_max - t_min)
        a_initial = log_y_min - b_initial * t_min
        return a_initial, b_initial

    def fit_slope(self, t, y):
        log_y = np.log(y)
        p0 = self.initial_guess(t, y)
        popt, _ = curve_fit(self.exp_decline, t, log_y, p0=p0)
        return popt

    def find_min_slope(self):
        unique_t = np.unique(self.t)
        if len(unique_t) < self.windowsize:
            raise ValueError(
                "Not enough unique time points for the specified window size."
            )

        min_slope = float("inf")
        best_fit_params = None
        best_t_fit = None
        best_y_fit = None

        for i in range(len(unique_t) - self.windowsize + 1):
            window_t = unique_t[i : i + self.windowsize]
            mask = np.isin(self.t, window_t)
            t_window = self.t[mask]
            y_window = self.y[mask]

            popt = self.fit_slope(t_window, y_window)
            slope = popt[1]

            if slope < min_slope:
                min_slope = slope
                best_fit_params = popt
                best_t_fit = t_window
                best_y_fit = y_window

        return best_fit_params, best_t_fit, best_y_fit

    def plot(self):
        plt.scatter(self.t, self.y, label="Data", color="blue")

        fitted_y = np.exp(self.popt[0] + self.popt[1] * self.t)
        plt.plot(self.t, fitted_y, label="Fitted Function", color="red")

        plt.xlabel("Time (t)")
        plt.ylabel("Signal (y)")
        plt.title("Exponential Decay Fit")
        plt.legend()
        plt.yscale("log")

    def summarize(self):
        return {
            "t": self.t,
            "y": self.y,
            "t_fit": self.t_fit,
            "y_fit": self.y_fit,
            "slope": self.popt[1],
            "intercept": self.popt[0],
        }
