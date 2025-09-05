# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome

"""

ppanggolin draw -p stepbystep/pangenome.h5 --tile_plot --nocloud --soft_core 0.92 --ucurve --output stepbystep -f
ppanggolin draw -p stepbystep/pangenome.h5 --draw_spots -o stepbystep -f
ppanggolin draw -p stepbystep/pangenome.h5 --draw_spots --spots all -o stepbystep -f
"""

OUTDIR_COMMANDS_AND_FILES = [
    (
        "draw_tile",
        "ppanggolin draw -p {pangenome} --tile_plot --nocloud --soft_core 0.92 "
        "--ucurve --output {outdir}",
        [
            "tile_plot.html",
            "Ushaped_plot.html",
        ],
    ),
    (
        "draw_spots",
        "ppanggolin draw -p {pangenome} --draw_spots -o {outdir} -f",
        [
            "spot_0.gexf",
            "spot_0.html",
            "spot_0_identical_rgps.tsv",
            "spot_1.gexf",
            "spot_1.html",
            "spot_1_identical_rgps.tsv",
        ],
    ),
    (
        "draw_spots_all",
        "ppanggolin draw -p {pangenome} --draw_spots --spots all -o {outdir} -f",
        [
            "spot_0.gexf",
            "spot_0.html",
            "spot_0_identical_rgps.tsv",
            "spot_1.gexf",
            "spot_1.html",
            "spot_1_identical_rgps.tsv",
        ],
    ),
]


@pytest.mark.parametrize(
    "test_outdir, cmd_template, expected_files",
    OUTDIR_COMMANDS_AND_FILES,
    ids=[case[0] for case in OUTDIR_COMMANDS_AND_FILES],
)
def test_stepbystep_outputs(
    stepbystep_pangenome, tmp_path, test_outdir, cmd_template, expected_files
):
    """Run each command on the prepared pangenome and check expected files exist."""
    outdir = tmp_path / test_outdir
    cmd = cmd_template.format(pangenome=stepbystep_pangenome, outdir=outdir)

    run_ppanggolin_command(cmd)

    for fname in expected_files:
        fpath = Path(outdir) / fname
        assert fpath.exists(), f"Expected file {fname} not found after `{cmd}`"
        assert fpath.stat().st_size > 0, f"File {fname} is empty after `{cmd}`"
