# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome

"""
ppanggolin write_genomes  -p stepbystep/pangenome.h5 --output stepbystep -f --fasta genomes.fasta.list --gff --proksee --table
"""

def test_write_genomes(stepbystep_pangenome):
    outdir = stepbystep_pangenome.parent / "write_genomes_outdir"
    pangenome = stepbystep_pangenome
    cmd = (
        f"ppanggolin write_genomes "
        f"-p {pangenome} "
        f"--output {outdir} "
        f"-f --fasta testingDataset/genomes.fasta.list "
        f"--gff --proksee --table"
    )

    run_ppanggolin_command(cmd)

    # Check the expected subdirectories exist
    expected_dirs = ["gff", "proksee", "table"]
    for sub in expected_dirs:
        subdir = outdir / sub
        assert subdir.exists() and subdir.is_dir(), f"Missing subdir {subdir}"

        # Check number of files inside
        files = [f for f in subdir.iterdir() if f.is_file()]
        assert len(files) == 51, f"Expected 51 files in {subdir}, found {len(files)}"
