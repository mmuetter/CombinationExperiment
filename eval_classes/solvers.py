solvers = {
    "Nelder-Mead": {
        "maxiter": 25000,  # Very high number to allow extensive searching
        "maxfev": 25000,  # Similarly high number for function evaluations
        "xatol": 1e-9,  # High precision for convergence in variable space
        "fatol": 1e-9,  # High precision for convergence in objective function value
    }
}
