from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_gbff_basic_wf import gbff_wf_pangenome
import logging

"""
ppanggolin msa --pangenome myannopang/pangenome.h5 --source dna --partition core -o myannopang/ -f --use_gene_id --phylo --single_copy --cpu $NUM_CPUS
"""


def test_msa(gbff_wf_pangenome, num_cpus, tmp_path):

    logging.info(gbff_wf_pangenome)

    msa_outdir = tmp_path / "msa_outdir"
    cmd = (
        f"ppanggolin msa --pangenome {gbff_wf_pangenome} --source dna "
        f"--partition core -o {msa_outdir} "
        f"--use_gene_id --phylo --single_copy --cpu {num_cpus}"
    )

    run_ppanggolin_command(cmd)

    subdir = msa_outdir / "msa_core_dna/"
    assert subdir.exists() and subdir.is_dir(), f"Missing subdir {subdir}"

    fname = "core_genome_alignment.aln"
    fpath = msa_outdir / fname
    assert fpath.exists(), f"Expected file {fname} not found after msa command: `{cmd}`"
    assert fpath.stat().st_size > 0, f"File {fname} is empty after msa command: `{cmd}`"
