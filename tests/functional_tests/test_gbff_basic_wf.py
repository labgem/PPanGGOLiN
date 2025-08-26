import pytest
import logging
from pathlib import Path

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
def gbff_wf_pangenome(tmp_path_factory, num_cpus, request):
    """Run ppanggolin workflow on GBFF files once and cache the result directory."""
    cache = request.config.cache
    cached_dir = cache.get("ppanggolin/gbff_wf_dir", None)

    if cached_dir and Path(cached_dir).exists():
        outdir = Path(cached_dir)
    else:
        tmp_path = tmp_path_factory.mktemp("pang_test")
        outdir = tmp_path / "myannopang"
        cmd = (
            f"ppanggolin workflow --cpu {num_cpus} "
            f"--anno testingDataset/genomes.gbff.list --output {outdir}"
        )
        run_ppanggolin_command(cmd)
        # Cache the directory for future pytest runs
        cache.set("ppanggolin/gbff_wf_dir", str(outdir))

    print(outdir)
    return outdir / "pangenome.h5"


#  cat myannopang/gene_families.tsv | cut -f1,2,4 > clusters.tsv
@pytest.fixture(scope="session")
def cluster_file(gbff_wf_pangenome, tmp_path_factory):

    gene_families_file = gbff_wf_pangenome.parent / "gene_families.tsv"

    cluster_file = tmp_path_factory.mktemp("cluster") / "clusters.tsv"

    with open(gene_families_file) as fin, open(cluster_file, "w") as fout:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                fout.write(f"{parts[0]}\t{parts[1]}\n")

    return cluster_file


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


def test_file_checksums(gbff_wf_pangenome, update_golden):
    pangenome_dir = gbff_wf_pangenome.parent
    file_path = pangenome_dir / "gene_families.tsv"
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)
