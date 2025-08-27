from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_fasta_all_wf import fasta_all_wf_pangenome
import logging

"""

ppanggolin align --pangenome mybasicpangenome/pangenome.h5 --sequences some_chlam_proteins.fasta \
                    --output test_align --draw_related --getinfo --fast --cpu $NUM_CPUS
                         
"""


def test_align(fasta_all_wf_pangenome, num_cpus, tmp_path):

    logging.info(fasta_all_wf_pangenome)

    align_outdir = tmp_path / "align_outdir"
    cmd = (
        f"ppanggolin align --pangenome {fasta_all_wf_pangenome} --sequences testingDataset/some_chlam_proteins.fasta"
        f" --output test_align --draw_related --getinfo --fast --cpu {num_cpus} --output {align_outdir}"
    )

    run_ppanggolin_command(cmd)

    expected_outputs = ["alignment_input_seqs_to_pangenome_gene_families.tsv",
                        "info_input_seq.tsv",
                        "sequences_partition_projection.tsv",
                        "spot_0.gexf",
                        "spot_0.html",
                        "spot_0_identical_rgps.tsv"]
    for fname in expected_outputs:
        fpath = align_outdir / fname
        assert fpath.exists(), f"Expected file {fname} not found after `{cmd}`"
        assert fpath.stat().st_size > 0, f"File {fname} is empty after `{cmd}`"