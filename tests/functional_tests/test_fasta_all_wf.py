import pytest

from pathlib import Path

from tests.utils.checksum import assert_or_update_hash
from tests.utils.file_compare import check_pangenome_info
from tests.utils.run_ppanggolin import run_ppanggolin_command
import logging
import shutil

# create a module-level logger
logger = logging.getLogger(__name__)


"""
ppanggolin all --cpu $NUM_CPUS --fasta genomes.fasta.list --output mybasicpangenome
ppanggolin info --pangenome mybasicpangenome/pangenome.h5 --content --parameters --status > info_to_test/mybasicpangenome_info.yaml
cat info_to_test/mybasicpangenome_info.yaml   
echo "$(grep 'mybasicpangenome/gene_families.tsv' expected_info_files/checksum.txt | cut -d' ' -f1)  mybasicpangenome/gene_families.tsv" | shasum -a 256 -c - || { echo 'Checksum verification failed.' >&2; exit 1; }
shasum -a 256 mybasicpangenome/gene_families.tsv > info_to_test/checksum.txt
cd -

"""


@pytest.fixture(scope="session")
def fasta_all_wf_pangenome(num_cpus, request):
    """Run ppanggolin all once and cache the result directory across pytest runs."""

    cache_dir = Path(request.config.cache.makedir("fasta_all_wf"))
    pangenome = cache_dir / "mybasicpangenome" / "pangenome.h5"

    if pangenome.exists():
        logger.warning(f"Reusing cached fasta_all_wf_pangenome: {pangenome}")
    else:
        logger.info(f"Recomputing fasta_all_wf_pangenome in {cache_dir}")
        # remove previous output directory if it exists
        outdir = cache_dir / "mybasicpangenome"
        if outdir.exists():
            shutil.rmtree(outdir)

        cmd = (
            f"ppanggolin all --cpu {num_cpus} "
            f"--fasta testingDataset/genomes.fasta.list --output {outdir}"
        )
        run_ppanggolin_command(cmd)

    return pangenome


def test_pangenome_created(fasta_all_wf_pangenome):
    assert fasta_all_wf_pangenome.exists()
    assert fasta_all_wf_pangenome.parent.exists()


def test_info_command(fasta_all_wf_pangenome, tmp_path, update_golden):
    info_file = tmp_path / "mybasicpangenome_info.yaml"
    check_pangenome_info(
        pangenome=fasta_all_wf_pangenome,
        info_file=info_file,
        update_golden=update_golden,
    )


def test_file_checksums(fasta_all_wf_pangenome, update_golden):
    pangenome_dir = fasta_all_wf_pangenome.parent
    file_path = pangenome_dir / "gene_families.tsv"
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)
