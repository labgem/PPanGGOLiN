import pytest
import logging

from tests.utils.checksum import assert_or_update_hash
from tests.utils.file_compare import check_pangenome_info
from tests.utils.run_ppanggolin import run_ppanggolin_command


"""

cd testingDataset
ppanggolin workflow --cpu $NUM_CPUS --anno genomes.gbff.list --output myannopang
ppanggolin info --pangenome myannopang/pangenome.h5 > info_to_test/myannopang_info.yaml
cat info_to_test/myannopang_info.yaml
echo "$(grep 'myannopang/gene_families.tsv' expected_info_files/checksum.txt | cut -d' ' -f1)  myannopang/gene_families.tsv" | shasum -a 256 -c - || { echo 'Checksum verification failed.' >&2; exit 1; }
shasum -a 256 myannopang/gene_families.tsv >> info_to_test/checksum.txt
"""

@pytest.fixture(scope="session")
def gbff_wf_pangenome(tmp_path_factory, num_cpus):
    tmp_path = tmp_path_factory.mktemp("pang_test")
    outdir = tmp_path / "myannopang"
    pangenome = outdir / "pangenome.h5"
    cmd = f"ppanggolin workflow --cpu {num_cpus} --anno testingDataset/genomes.gbff.list  --output {outdir}"

    run_ppanggolin_command(cmd)

    return pangenome


def test_pangenome_created(gbff_wf_pangenome):
    
    logging.info(gbff_wf_pangenome)
    

    pangenome_dir = gbff_wf_pangenome.parent
    assert pangenome_dir.exists()
    assert gbff_wf_pangenome.exists()


def test_info_command(gbff_wf_pangenome, tmp_path, update_golden):
    info_file = tmp_path / "myannopang_info.yaml"
    check_pangenome_info(
        pangenome=gbff_wf_pangenome, info_file=info_file, update_golden=update_golden
    )


FILES_TO_CHECK = ["gene_families.tsv"]
@pytest.mark.parametrize("fname", FILES_TO_CHECK)
def test_file_checksums(gbff_wf_pangenome, fname, update_golden):
    pangenome_dir = gbff_wf_pangenome.parent
    file_path = pangenome_dir / fname
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)


