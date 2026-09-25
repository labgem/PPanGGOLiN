# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import networkx as nx
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome

OUTDIR_COMMANDS_AND_FILES = [
    (
        "write_pangenome_outdir",
        "ppanggolin write_pangenome -p {pangenome} --output {outdir} -f "
        "--soft_core 0.9 --dup_margin 0.06 --gexf --light_gexf --gt --csv --Rtab --stats "
        "--partitions --compress --json --spots --regions --borders --families_tsv --cpu 1",
        [
            "border_protein_genes.fasta.gz",
            "genomes_statistics.tsv.gz",
            "pangenomeGraph.gexf.gz",
            "pangenomeGraph.gt.gz",
            "gene_families.tsv.gz",
            "matrix.csv.gz",
            "pangenomeGraph.json.gz",
            "regions_of_genomic_plasticity.tsv.gz",
            "summarize_spots.tsv.gz",
            "gene_presence_absence.Rtab.gz",
            "mean_persistent_duplication.tsv.gz",
            "pangenomeGraph_light.gexf.gz",
            "spot_borders.tsv.gz",
            "spots.tsv.gz",
            "partitions/cloud.txt",
            "partitions/exact_accessory.txt",
            "partitions/exact_core.txt",
            "partitions/persistent.txt",
            "partitions/S1.txt",
            "partitions/shell.txt",
            "partitions/soft_accessory.txt",
            "partitions/soft_core.txt",
        ],
    ),
    (
        "write_pangenome_outdir_json",
        "ppanggolin write_pangenome -p {pangenome} --output {outdir} --json",
        [
            "pangenomeGraph.json",
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

        if ".gexf" in fname:
            # Existence is not enough: a malformed header still writes a
            # non-empty file that no namespace-aware parser can read.
            try:
                graph = nx.read_gexf(fpath)
            except Exception as err:
                pytest.fail(f"{fname} was written but could not be parsed back: {err}")
            assert graph.number_of_nodes() > 0, f"{fname} parsed but has no nodes"
            # Fixing only the default namespace and leaving xmlns:viz behind
            # still parses, but drops every colour and size.
            assert any(
                "viz" in data for _, data in graph.nodes(data=True)
            ), f"{fname} parsed but its viz attributes were dropped"
