from pypetting_extra import default_worktable as worktable
from labware import labwares
import numpy as np

carrier = worktable.carrier

storex = worktable.incubator
shelf = worktable.carrier["Shelf 8x4Pos"]

tips_I = worktable.carrier["MCA96 3Pos"].define_labware(labwares["diti"], 1)


mca = worktable.mca
mca.tips = [tips_I]
tip_arr1 = np.array([False, True, True, True] + 4 * [False])
tip_arr2 = np.array(4 * [False] + [True, True, True, False])
tip_arr = tip_arr1 + tip_arr2


liha = worktable.liha


roma = worktable.roma
tilter = worktable.tilter
lid1_pos = worktable.carrier["MP 2Pos Fixed"].gridsite(0)
lid2_pos = worktable.carrier["MP 2Pos Fixed"].gridsite(1)
lid3_pos = shelf.gridsite(30)
rotated_site = tilter.gridsite(0)
plate1_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(0)
plate2_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(1)
plate3_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(2)
shaker_pos = worktable.carrier["Te-Shake"].gridsite(0)
