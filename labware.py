from pypetting import Labware
from pypetting.labware import labwares

corning_8row = Labware("8 Row DeepWell Corning", 8, 1)
deepwell96 = Labware("96 DeepWell Greiner", 8, 12)
dw12col = Labware("12Col Through", 8, 12)
rotated_12col = Labware("12Col rotated", 12, 8)
greiner96 = labwares["greiner96"]
greiner384 = labwares["greiner384"]
trough300 = Labware("Trough 300ml MCA96", 8, 12)


labwares.update(
    {
        "corning_8row": corning_8row,
        "deepwell96": deepwell96,
        "dw12col": dw12col,
        "rotated_12col": rotated_12col,
        "diti": Labware("DiTi 200ul SBS MCA96", 8, 12),
        "trough300": trough300,
    }
)
