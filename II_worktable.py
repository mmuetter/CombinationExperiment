from pypetting_extra import default_worktable as worktable
from labware import labwares
import numpy as np


storex = worktable.incubator
shelf = worktable.carrier["Shelf 8x4Pos"]

tips_I = worktable.carrier["MP 3Pos Deck"].define_labware(labwares["diti"], 0)
tips_II = worktable.carrier["MP 3Pos Deck"].define_labware(labwares["diti"], 1)

worktable.add_tips(tips_I)
worktable.add_tips(tips_II)
tip_arr1 = np.array([False, True, True, True] + 4 * [False])
tip_arr2 = np.array(4 * [False] + [True, True, True, False])
tip_arr = tip_arr1 + tip_arr2

liha = worktable.liha

mca = worktable.mca
roma = worktable.roma
tilter = worktable.tilter
lid1_pos = worktable.carrier["MP 2Pos Fixed"].gridsite(0)
lid2_pos = worktable.carrier["MP 2Pos Fixed"].gridsite(1)
lid3_pos = shelf.gridsite(30)
rotated_site = tilter.gridsite(0)
plate1_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(0)
plate2_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(1)
plate3_pos = worktable.carrier["MP 3Pos Fixed"].gridsite(2)
