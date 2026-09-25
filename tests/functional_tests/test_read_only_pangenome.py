import os
import shutil

import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome


def test_read_only_pangenome(stepbystep_pangenome, tmp_path):
    """Commands that only read a pangenome must not require write access to it.

    Reference collections are distributed on read-only media, and making a
    local copy read-only is the simplest way to keep a pangenome byte-identical
    to the one it was downloaded from.
    """
    pangenome = tmp_path / "pangenome.h5"
    shutil.copy(stepbystep_pangenome, pangenome)
    pangenome.chmod(0o444)

    if os.access(pangenome, os.W_OK):
        # Running as root, or on a filesystem that ignores the mode. PyTables
        # gates write modes on this very call, so the commands below would
        # succeed even unfixed and the test would prove nothing.
        pytest.skip("write access to the pangenome could not be revoked")

    run_ppanggolin_command(f"ppanggolin info -p {pangenome}")

    # Genome fluidity is already stored in the file, so this reads it back
    # and prints it without recomputing anything.
    run_ppanggolin_command(f"ppanggolin metrics -p {pangenome} --genome_fluidity")
