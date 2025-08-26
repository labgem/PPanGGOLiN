# tests/functional_tests/test_stepbystep.py
from logging import config
from pathlib import Path

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome

from tests.functional_tests.test_fasta_all_wf import fasta_all_wf_pangenome
from tests.functional_tests.test_gbff_basic_wf import cluster_file, gbff_wf_pangenome

from tests.utils.checksum import assert_or_update_hash
from tests.utils.file_compare import check_pangenome_info
from tests.utils.checksum import assert_or_update_hash


import logging

import pytest
"""
ppanggolin utils --default_config panrgp -o panrgp_default_config.yaml
cut -f1,2 clusters.tsv > clusters_without_frag.tsv
ppanggolin panrgp  --anno genomes.gbff.list --cluster clusters_without_frag.tsv -o test_config --config panrgp_default_config.yaml --cpu $NUM_CPUS
echo "$(grep 'test_config/gene_families.tsv' expected_info_files/checksum.txt | cut -d' ' -f1)  test_config/gene_families.tsv" | shasum -a 256 -c - || { echo 'Checksum verification failed.' >&2; exit 1; }
shasum -a 256 test_config/gene_families.tsv >> info_to_test/checksum.txt
"""


@pytest.fixture(scope="session")
def pangenome_default_config(cluster_file, tmp_path_factory, num_cpus, request):
    """Run ppanggolin all once and cache the result across pytest runs."""
    cache = request.config.cache
    cached_path = cache.get("ppanggolin/config_file_pangenome", None)

    if cached_path and Path(cached_path).exists():
        # reuse cached result
        return Path(cached_path)

    # Otherwise recompute
    tmp_path = tmp_path_factory.mktemp("pang_config_test")
    outdir = tmp_path / "config_pangenome"
    pangenome = outdir / "pangenome.h5"
    config_file = tmp_path / "panrgp_default_config.yaml"
    
    cmd = f"ppanggolin utils --default_config panrgp -o {config_file}"
    run_ppanggolin_command(cmd)

    
    cmd = (
        f"ppanggolin panrgp --cpu {num_cpus} "
        f"--anno testingDataset/genomes.gbff.list --output {outdir} "
        f"--cluster {cluster_file} -o {outdir} --config {config_file} "
        f"--cpu {num_cpus}"
    )
    
    run_ppanggolin_command(cmd)

    # Save to pytest cache for reuse in later pytest runs
    cache.set("ppanggolin/config_file_pangenome", str(pangenome))
    cache.set("ppanggolin/config_file_pangenome_gene_families.tsv", str(outdir / "gene_families.tsv"))

    return pangenome




def test_pangenome_created(pangenome_default_config):
    
    logging.info(pangenome_default_config)
    

    pangenome_dir = pangenome_default_config.parent
    assert pangenome_dir.exists()
    assert pangenome_default_config.exists()


FILES_TO_CHECK = ["gene_families.tsv"]
@pytest.mark.parametrize("fname", FILES_TO_CHECK)
def test_file_checksums(pangenome_default_config, fname, update_golden):
    pangenome_dir = pangenome_default_config.parent
    file_path = pangenome_dir / fname
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)
