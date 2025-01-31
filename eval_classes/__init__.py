from .bootstrap import BootStrapRamp
from .functions import (
    call,
    collect,
    violine_plot,
    t_plot,
    compare_methods,
    get_permutation,
    get_rep_rates,
    plot_rep_rates,
    permute_rates,
    combine_signals,
)
from .permutation import PermutationTest

# from .experiment import Experiment, PathManager
from .sec_selection import SelectSection
from .solvers import solvers
from .plots import cfu_rlu_timelplot
from .import_lum import ImportLumMeasures
from .timeplot import plot, timeplot
from .full_eval import full_evaluation
from .recover import recover_sec_df
from .paperfig import plot_df
from .bootstrap2 import BootstrapSlope
from .compare_pi import CompareSlopes
from .window_fit import WindowFit
from .fit_slope import FitSlope
