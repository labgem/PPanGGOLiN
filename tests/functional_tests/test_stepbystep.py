# tests/functional_tests/test_stepbystep.py
from pathlib import Path
from unittest.mock import patch
import pytest

from tests.utils.checksum import assert_or_update_hash
from tests.utils.file_compare import check_pangenome_info
from tests.utils.run_ppanggolin import run_ppanggolin_command


@pytest.fixture(scope="module")
def outdir(tmp_path_factory):
    return tmp_path_factory.mktemp("stepbystep_test") / "stepbystep"


@pytest.fixture(scope="module")
def annotate(outdir, num_cpus):
    cmd = f"ppanggolin annotate --fasta testingDataset/genomes.fasta.list --output {outdir} --kingdom bacteria --cpu {num_cpus}"
    run_ppanggolin_command(cmd)
    return outdir


def test_annotate(annotate):
    pangenome = annotate / "pangenome.h5"
    assert pangenome.exists()


@pytest.fixture(scope="module")
def cluster(annotate, num_cpus):
    pangenome = annotate / "pangenome.h5"
    cmd = f"ppanggolin cluster -p {pangenome} --coverage 0.8 --identity 0.8 --cpu {num_cpus}"
    run_ppanggolin_command(cmd)
    return annotate


@pytest.fixture(scope="module")
def graph(cluster):
    pangenome = cluster / "pangenome.h5"
    cmd = f"ppanggolin graph -p {pangenome} -r 10"
    run_ppanggolin_command(cmd)
    return cluster


@pytest.fixture(scope="module")
def partition(graph, num_cpus):
    pangenome = graph / "pangenome.h5"
    outdir = graph
    cmd = f"ppanggolin partition --output {outdir} -f -p {pangenome} --cpu {num_cpus} -b 2.6 -ms 10 -fd -ck 500 -Kmm 3 12 -im 0.04 --draw_ICL"
    print(cmd)
    run_ppanggolin_command(cmd)
    return outdir


@pytest.fixture(scope="module")
def rgp(partition):
    pangenome = partition / "pangenome.h5"
    outdir = partition
    cmd = f"ppanggolin rgp -p {pangenome} --persistent_penalty 2 --variable_gain 1 --min_score 3 --dup_margin 0.05"
    run_ppanggolin_command(cmd)
    return outdir


@pytest.fixture(scope="module")
def spot(rgp):
    pangenome = rgp / "pangenome.h5"
    outdir = rgp
    cmd = f"ppanggolin spot -p {pangenome} --output {outdir} --spot_graph --overlapping_match 2 --set_size 3 --exact_match_size 1 -f"
    run_ppanggolin_command(cmd)
    return outdir


@pytest.fixture(scope="module")
def module(spot):
    pangenome = spot / "pangenome.h5"
    outdir = spot
    cmd = f"ppanggolin module -p {pangenome} --transitive 4 --size 3 --jaccard 0.86 --dup_margin 0.05"
    run_ppanggolin_command(cmd)
    return outdir

@pytest.fixture(scope="module")
def metrics(module, num_cpus):
    cmd = f"ppanggolin metrics -p {module / 'pangenome.h5'} --genome_fluidity --no_print_info --recompute_metrics --log metrics.log"
    run_ppanggolin_command(cmd)
    return module


def test_cluster(cluster):
    pangenome = cluster / "pangenome.h5"
    assert pangenome.exists()


def test_graph(graph):
    pangenome = graph / "pangenome.h5"
    assert pangenome.exists()


def test_partition(partition):
    # check partition outputs exist
    pangenome = partition / "pangenome.h5"
    assert pangenome.exists()


def test_rgp(rgp):
    pangenome = rgp / "pangenome.h5"
    assert pangenome.exists()


def test_spot(spot):
    pangenome = spot / "pangenome.h5"
    assert pangenome.exists()


def test_module(module):
    pangenome = module / "pangenome.h5"
    assert pangenome.exists()


def test_rarefaction(module, num_cpus):
    outdir = module
    cmds = [
        f"ppanggolin rarefaction --output {outdir} -f -p {outdir / 'pangenome.h5'} --depth 5 --min 1 --max 50 -ms 10 -fd -ck 30 -K 3 --soft_core 0.9 -se 42 --cpu {num_cpus}",
        f"ppanggolin rarefaction --output {outdir} -f -p {outdir / 'pangenome.h5'} --reestimate_K --depth 1 --min 40 --cpu {num_cpus}",
    ]
    for cmd in cmds:
        run_ppanggolin_command(cmd)


def test_draw_and_write(module, update_golden):
    outdir = module
    cmds = [
        f"ppanggolin draw -p {outdir / 'pangenome.h5'} --tile_plot --nocloud --soft_core 0.92 --ucurve --output {outdir} -f",
        f"ppanggolin draw -p {outdir / 'pangenome.h5'} --draw_spots -o {outdir} -f",
        f"ppanggolin write_pangenome -p {outdir / 'pangenome.h5'} --output {outdir} -f --soft_core 0.9 --dup_margin 0.06 --gexf --light_gexf --csv --Rtab --stats --partitions --compress --json --spots --regions --borders --families_tsv --cpu 1",
        f"ppanggolin write_genomes -p {outdir / 'pangenome.h5'} --output {outdir} -f --fasta testingDataset/genomes.fasta.list --gff --proksee --table",
    ]
    for cmd in cmds:
        run_ppanggolin_command(cmd)

    file_path = outdir / "gene_families.tsv.gz"
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)


def test_info_command(metrics, tmp_path, update_golden):
    info_file = tmp_path / "stepbystep_info.yaml"
    pangenome = metrics / "pangenome.h5"
    check_pangenome_info(
        pangenome=pangenome, info_file=info_file, update_golden=update_golden
    )
