from tests.functional_tests.test_fasta_all_wf import fasta_all_wf_pangenome
from tests.functional_tests.test_metadata_pan import pangenome_with_metadata

# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command


"""

ppanggolin rgp_cluster --pangenome mybasicpangenome/pangenome.h5
ppanggolin rgp_cluster --pangenome mybasicpangenome/pangenome.h5 --ignore_incomplete_rgp --grr_metric max_grr -f --graph_formats graphml gexf
ppanggolin rgp_cluster --pangenome mybasicpangenome/pangenome.h5 --no_identical_rgp_merging -o rgp_clustering_no_identical_rgp_merging --graph_formats graphml

"""


OUTDIR_COMMANDS_AND_FILES = [
    (
        "rgp_cluster_default",
        "ppanggolin rgp_cluster --pangenome {pangenome} --output {outdir}",
        ["rgp_cluster.gexf", "rgp_cluster.tsv"],
    ),
    (
        "rgp_cluster_graph_formats",
        "ppanggolin rgp_cluster --pangenome {pangenome} "
        "--ignore_incomplete_rgp --grr_metric max_grr -f "
        "--graph_formats graphml gexf --output {outdir}",
        ["rgp_cluster.gexf", "rgp_cluster.graphml", "rgp_cluster.tsv"],
    ),
    (
        "rgp_clustering_no_identical_rgp_merging",
        "ppanggolin rgp_cluster --pangenome {pangenome} "
        "--no_identical_rgp_merging --graph_formats graphml "
        "--output {outdir}",
        ["rgp_cluster.graphml", "rgp_cluster.tsv"],
    ),
]


@pytest.mark.parametrize(
    "test_outdir, cmd_template, expected_files",
    OUTDIR_COMMANDS_AND_FILES,
    ids=[case[0] for case in OUTDIR_COMMANDS_AND_FILES],
)
def test_rgp_cluster_outputs(
    fasta_all_wf_pangenome, tmp_path, test_outdir, cmd_template, expected_files
):
    """Run rgp_cluster with various options and check expected files exist."""
    outdir = tmp_path / test_outdir
    cmd = cmd_template.format(pangenome=fasta_all_wf_pangenome, outdir=outdir)

    run_ppanggolin_command(cmd)

    for fname in expected_files:
        fpath = Path(outdir) / fname
        assert fpath.exists(), f"Expected file {fname} not found after `{cmd}`"
        assert fpath.stat().st_size > 0, f"File {fname} is empty after `{cmd}`"


def test_rgp_cluster_with_metadata(pangenome_with_metadata, tmp_path):
    """Run rgp_cluster with various options and check expected files exist."""

    outdir = tmp_path / "rgp_cluster_with_metadata"

    cmd = (
        f"ppanggolin rgp_cluster --pangenome {pangenome_with_metadata} "
        f"--output {outdir} --add_metadata"
    )
    run_ppanggolin_command(cmd)

    expected_files = ["rgp_cluster.gexf", "rgp_cluster.tsv", "rgp_cluster.graphml"]
    for fname in expected_files:
        fpath = Path(outdir) / fname
        assert fpath.exists(), f"Expected file {fname} not found after `{cmd}`"
        assert fpath.stat().st_size > 0, f"File {fname} is empty after `{cmd}`"
