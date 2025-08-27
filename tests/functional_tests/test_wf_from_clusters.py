import pytest
import logging
from pathlib import Path

from tests.utils.checksum import assert_or_update_hash
from tests.utils.file_compare import check_pangenome_info
from tests.utils.run_ppanggolin import run_ppanggolin_command

from tests.functional_tests.test_gbff_basic_wf import (
    cluster_file,
    cluster_file_with_fragments,
    cluster_file_with_representative,
    gbff_wf_pangenome,
)
from tests.utils.cache_utils import run_with_cache


logger = logging.getLogger(__name__)

"""
cd testingDataset
cat myannopang/gene_families.tsv | cut -f1,2,4 > clusters.tsv 
ppanggolin panrgp --anno genomes.gbff.list --cluster clusters.tsv --output readclusterpang  --cpu $NUM_CPUS 
ppanggolin annotate --anno genomes.gbff.list --output readclusters --cpu $NUM_CPUS
awk 'BEGIN{FS=OFS="\t"} {$1 = $1 OFS $1} 1' clusters.tsv > clusters_with_reprez.tsv;
ppanggolin cluster --clusters clusters_with_reprez.tsv -p readclusters/pangenome.h5 --cpu $NUM_CPUS
ppanggolin msa --pangenome readclusterpang/pangenome.h5 --partition persistent --phylo -o readclusterpang/msa/ -f --cpu $NUM_CPUS
echo "$(grep 'readclusterpang/gene_families.tsv' expected_info_files/checksum.txt | cut -d' ' -f1)  readclusterpang/gene_families.tsv" | shasum -a 256 -c - || { echo 'Checksum verification failed.' >&2; exit 1; }
shasum -a 256 readclusterpang/gene_families.tsv >> info_to_test/checksum.txt
cd -
"""


@pytest.fixture(scope="session")
def gbff_panrgp_from_cluster_pangenome(cluster_file, num_cpus, request):

    cache_dir = Path(request.config.cache.makedir("pang_test"))
    outdir = cache_dir / "readclusterpang"

    cmd = (
        f"ppanggolin panrgp --cpu {num_cpus} "
        f"--anno testingDataset/genomes.gbff.list --output {outdir} "
        f"--cluster {cluster_file}"
    )

    outdir = run_with_cache(
        request,
        cache_key="ppanggolin/gbff_wf_dir",
        outdir=outdir,
        cmds=[cmd],
    )

    return outdir / "pangenome.h5"


def test_pangenome_created(gbff_panrgp_from_cluster_pangenome):

    logging.info(gbff_panrgp_from_cluster_pangenome)

    pangenome_dir = gbff_panrgp_from_cluster_pangenome.parent
    assert pangenome_dir.exists()
    assert gbff_panrgp_from_cluster_pangenome.exists()


def test_info_command(gbff_panrgp_from_cluster_pangenome, tmp_path, update_golden):
    info_file = tmp_path / "panrgp_from_existing_cluster_info.yaml"
    check_pangenome_info(
        pangenome=gbff_panrgp_from_cluster_pangenome,
        info_file=info_file,
        update_golden=update_golden,
    )


def test_file_checksums(gbff_panrgp_from_cluster_pangenome, update_golden):
    pangenome_dir = gbff_panrgp_from_cluster_pangenome.parent
    file_path = pangenome_dir / "gene_families.tsv"
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)


@pytest.fixture(scope="session")
def annotate_pangenome(tmp_path_factory, num_cpus):

    tmp_path = tmp_path_factory.mktemp("pang_test")
    outdir = tmp_path / "readclusters"
    cmd = (
        f"ppanggolin annotate --cpu {num_cpus} "
        f"--anno testingDataset/genomes.gbff.list --output {outdir} "
    )
    run_ppanggolin_command(cmd)

    return outdir / "pangenome.h5"


def test_cluster_from_existing_clustering(
    annotate_pangenome, cluster_file_with_representative, update_golden
):

    cmd = f"ppanggolin cluster  -p {annotate_pangenome} --clusters {cluster_file_with_representative}"
    outdir = annotate_pangenome.parent
    run_ppanggolin_command(cmd)

    cmd = f"ppanggolin write_pangenome -p {annotate_pangenome} --output {outdir} -f --families_tsv -f"

    run_ppanggolin_command(cmd)

    file_path = outdir / "gene_families.tsv"
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)
