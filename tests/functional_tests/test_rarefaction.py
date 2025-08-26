from pathlib import Path
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome

"""
ppanggolin rarefaction --output stepbystep -f -p stepbystep/pangenome.h5 --depth 5 --min 1 --max 50 -ms 10 -fd -ck 30 -K 3 --soft_core 0.9 -se $RANDOM

# Check reestimate K with depth 1 and min 40 to avoid long computation time
ppanggolin rarefaction --output stepbystep -f -p stepbystep/pangenome.h5 --reestimate_K --depth 1 --min 40

"""

OUTDIR_COMMANDS_AND_FILES = [
    (
        "rarefaction1",
        "ppanggolin rarefaction --output {outdir} -f -p {pangenome} "
        "--depth 2 --min 1 --max 50 -ms 10 -fd -ck 30 -K 3 --soft_core 0.9 -se 42",
        [
            "rarefaction.csv",
            "rarefaction_curve.html",
            "rarefaction_parameters.csv",
        ],
    ),
    (
        "rarefaction_reestimate_K",
        "ppanggolin rarefaction --output {outdir} -f -p {pangenome} "
        "--reestimate_K --depth 1 --min 40",
        [
            "rarefaction.csv",
            "rarefaction_curve.html",
            "rarefaction_parameters.csv",
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
