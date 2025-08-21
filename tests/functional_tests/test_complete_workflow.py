import pytest


from tests.utils.checksum import assert_or_update_hash
from tests.utils.file_compare import check_pangenome_info
from tests.utils.run_ppanggolin import run_ppanggolin_command


@pytest.fixture(scope="module")
def pangenome_dir(tmp_path_factory, num_cpus):
    tmp_path = tmp_path_factory.mktemp("pang_test")
    outdir = tmp_path / "mybasicpangenome"

    cmd = f"ppanggolin all --cpu {num_cpus} --fasta testingDataset/genomes.fasta.list --output {outdir}"
    run_ppanggolin_command(cmd)

    return outdir


def test_pangenome_created(pangenome_dir):
    assert pangenome_dir.exists()
    assert (pangenome_dir / "pangenome.h5").exists()


def test_info_command(pangenome_dir, tmp_path, update_golden):
    info_file = tmp_path / "mybasicpangenome_info.yaml"
    pangenome = pangenome_dir / "pangenome.h5"
    check_pangenome_info(
        pangenome=pangenome, info_file=info_file, update_golden=update_golden
    )


FILES_TO_CHECK = ["gene_families.tsv"]


@pytest.mark.parametrize("fname", FILES_TO_CHECK)
def test_file_checksums(pangenome_dir, fname, update_golden):
    file_path = pangenome_dir / fname
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)
