from pypetting_extra import Experiment
from configurations import AssayConfiguration
from workflow import AssayWorkflow
from plate_files import PlateFiles
from pypetting import GridSite, Labware, aspirate
from pypetting.labware import labwares
from setup import AssaySetup, storex, worktable, plate1_pos, plate2_pos, liha, roma


folder = "validation_experiments_2024"
exp_name = "test4"
exp_path = "/Users/malte/polybox/Shared/Robot-Malte/Luminescence/" + folder
experiment = Experiment(
    exp_name,
    exp_path,
    "C:\\Users\\COMPUTER\\polybox\\Robot-Malte\\Luminescence\\" + folder,
)

config = AssayConfiguration()
setup = AssaySetup(config, experiment)


setup.define_antibiotic_plate()
setup.define_plates()
setup.define_overnight_plate()
setup.define_PBS_trough()


def move_liha(
    grid_site: GridSite,
    column: int,
    positions: ArrayLike = None,
    spacing: int = 1,
    labware: Labware | str = "greiner96",
    local: bool = False,
    z_pos: int = 0,
    speed: int = 10,
):
    """Move LiHa to grid site."""

    if isinstance(labware, str):
        labware = labwares[labware]

    if positions is None:
        positions = np.array([True] * 8)

    if column not in range(1, labware.cols + 1):
        raise IndexError(f"Plate column out of bounds: {column=}")

    command = (
        (
            "B;MoveLiha("
            "255,"
            f"{grid_site.grid},"
            f"{grid_site.site},"
            f'{spacing},"'
        ).encode()
        + _well_select(positions, column, labware.rows, labware.cols)
        + (f'",{local:#d},{z_pos},0,{speed},0,0);').encode()
    )

    return command


plate = worktable.carrier["MP 3Pos Fixed"].define_plate(
    "test", labwares["greiner384"], 1
)


wl = experiment.setup_worklist("test.gwl")
column_mask = 6 * [False, True] + 2 * [False, True]

wl.add(liha.aspirate(plate, 1, 10, column_mask, retract=False))
wl.add([move_liha(plate.gridsite, column=1, positions=column_mask, local=True)])
wl.add(roma.move_plate(plate, plate1_pos))
