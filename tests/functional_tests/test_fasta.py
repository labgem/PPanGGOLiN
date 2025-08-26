# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome

"""
ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --prot_families all --gene_families shell --regions all --fasta genomes.fasta.list
ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --prot_families rgp --gene_families rgp --compress 
ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --prot_families softcore --gene_families softcore 
ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --prot_families module_0
ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --genes core --proteins cloud
ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --gene_families module_0 --genes module_0 --compress
ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --proteins cloud --cpu $NUM_CPUS --keep_tmp --compress

"""

OUTDIR_COMMANDS_AND_FILES = [
    (
        "fasta_rgp",
        "ppanggolin fasta -p {pangenome} --output {outdir} -f "
        "--prot_families rgp --gene_families rgp --compress",
        [
            "rgp_nucleotide_families.fasta.gz",
            "rgp_protein_families.faa.gz",
        ],
    ),
    (
        "fasta_softcore",
        "ppanggolin fasta -p {pangenome} --output {outdir} -f "
        "--prot_families softcore --gene_families softcore",
        [
            "softcore_nucleotide_families.fasta",
            "softcore_protein_families.faa",
        ],
    ),
    (
        "fasta_module0",
        "ppanggolin fasta -p {pangenome} --output {outdir} -f --prot_families module_0",
        [
            "module_0_protein_families.faa",
        ],
    ),
    (
        "fasta_core_cloud",
        "ppanggolin fasta -p {pangenome} --output {outdir} -f --genes core --proteins cloud",
        [
            "cloud_protein_genes.faa",
            "core_genes.fna",
        ],
    ),
    (
        "fasta_module0_compress",
        "ppanggolin fasta -p {pangenome} --output {outdir} -f "
        "--gene_families module_0 --genes module_0 --compress",
        [
            "module_0_genes.fna.gz",
            "module_0_nucleotide_families.fasta.gz",
        ],
    ),
    (
        "fasta_cloud_compress",
        "ppanggolin fasta -p {pangenome} --output {outdir} -f "
        "--proteins cloud --cpu 2 --keep_tmp --compress",
        [
            "cloud_protein_genes.faa.gz",
        ],
    ),
    (
        "fasta_all_shell",
        "ppanggolin fasta -p {pangenome} --output {outdir} -f "
        "--prot_families all --gene_families shell --regions all "
        "--fasta testingDataset/genomes.fasta.list",
        [
            "all_protein_families.faa",
            "all_rgp_genomic_sequences.fasta",
            "shell_nucleotide_families.fasta",
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
