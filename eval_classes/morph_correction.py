import numpy as np
from scipy.optimize import curve_fit


def linear_func(x, a, b):
    return a * x + b


def fit_v_params(morph_summary, antibiotic, v_pre):
    t = np.array([0, 2])
    v_past = morph_summary.loc[antibiotic, "volume"]
    y = np.array([v_pre, v_past])
    params, _ = curve_fit(linear_func, t, y)
    slope, intercept = params
    return slope, intercept


def morph_correction(t, y, antibiotic, morph_summary):
    v_pre = morph_summary.loc["control", "volume"]
    params = fit_v_params(morph_summary, antibiotic, v_pre)
    v_t = linear_func(t, *params)
    if params[0] > 0:
        return v_pre / v_t * y
    else:
        return y
