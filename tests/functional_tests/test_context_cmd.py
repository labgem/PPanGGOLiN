# tests/functional_tests/test_stepbystep.py
from pathlib import Path
from venv import logger
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_fasta_all_wf import fasta_all_wf_pangenome
from tests.functional_tests.test_gbff_basic_wf import gbff_wf_pangenome
import logging

logger = logging.getLogger(__name__)

"""
"""


def test_context_by_sequences(gbff_wf_pangenome, num_cpus, tmp_path):
    """Test ppanggolin context using input sequences."""

    logger.info("Running context test with sequences")
    
    outdir = tmp_path / "test_context"
    cmd = (
        f"ppanggolin context --pangenome {gbff_wf_pangenome} "
        f"--sequences testingDataset/some_chlam_proteins.fasta "
        f"--output {outdir} --fast --cpu {num_cpus}"
    )

    run_ppanggolin_command(cmd)

    expected_files = [
        "alignment_input_seqs_to_pangenome_gene_families.tsv",
        "gene_contexts.tsv",
        "gene_to_gene_family.tsv",
        "graph_context.graphml",
        "sequences_partition_projection.tsv"
    ]

    for fname in expected_files:
        fpath = outdir / fname
        assert fpath.exists(), f"Expected file {fname} not found after `{cmd}`"
        assert fpath.stat().st_size > 0, f"File {fname} is empty after `{cmd}`"


def test_context_by_family_id(gbff_wf_pangenome, num_cpus, tmp_path):
    """Test ppanggolin context using a single gene family ID."""

    logger.info("Running context test with a gene family ID")

    # Prepare input file with one family
    family_file = tmp_path / "one_family_of_module_1.txt"
    family_file.write_text("AP288_RS05055\n")

    outdir = tmp_path / "test_context_from_id"
    cmd = (
        f"ppanggolin context --pangenome {gbff_wf_pangenome} "
        f"--family {family_file} --output {outdir} --cpu {num_cpus}"
    )

    run_ppanggolin_command(cmd)

    expected_files = [
        "gene_contexts.tsv",
        "graph_context.graphml" 
    ]

    for fname in expected_files:
        fpath = outdir / fname
        assert fpath.exists(), f"Expected file {fname} not found after `{cmd}`"
        assert fpath.stat().st_size > 0, f"File {fname} is empty after `{cmd}`"
