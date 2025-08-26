# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome

OUTDIR_COMMANDS_AND_FILES = [
    (
        "write_pangenome_outdir",
        "ppanggolin write_pangenome -p {pangenome} --output {outdir} -f "
        "--soft_core 0.9 --dup_margin 0.06 --gexf --light_gexf --csv --Rtab --stats "
        "--partitions --compress --json --spots --regions --borders --families_tsv --cpu 1",
        [
            "border_protein_genes.fasta.gz",
            "genomes_statistics.tsv.gz",
            "pangenomeGraph.gexf.gz",
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
