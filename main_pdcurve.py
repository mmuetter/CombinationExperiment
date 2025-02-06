from pypetting_extra import Experiment
from configurations import AssayConfiguration
from setup_pdcurve import pdSetup
from workflow_pdcurve import pdWorkflow
from general_classes import PathManager
import pandas as pd
import numpy as np

## use pyenv 3.12.0

folder = "twofold1to1"
exp_name = "0_025_32_06022025"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/CombinationProject/" + folder

## The order of the plates array defines the order in which the plates are created.
# This is important, low concentrations should be first, because the same tips are used.
plates = ["antibiotic0", "antibiotic025", "antibiotic32"]
# the read order defines the order in which the plates are read.
# Plates may appear multiple times in the read order, but only once in the plates array.
read_order = ["assay32", "assay0", "assay32", "assay025"]


pm = PathManager(basepath="current")
drugs = pd.read_excel(pm.file_path("drugs.xlsx", folder="notes"))

experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\CombinationProject\\" + folder,
)

config = AssayConfiguration()
setup = pdSetup(config, experiment, plates)

workflow = pdWorkflow(setup)
workflow.grow_overnight()
workflow.prepare_exp_cultures()
workflow.fill_replicates()
workflow.treat_cultures()
workflow.luminescence_read_loop(read_order)

## EVOSCRIPT
# WL-Grow Overnight
# WL-Prepare exponential cultures
# WL-Treat cultures
# for i:
#     WL-Lumread plate i

wl = workflow.setup_worklist("lum_read_test.gwl")
wl.add(setup.lumread.measure("test.xml"))
