# tests/functional_tests/test_stepbystep.py
from pathlib import Path
from venv import logger
import pytest
import logging

from tests.utils.checksum import assert_or_update_hash
from tests.utils.file_compare import check_pangenome_info
from tests.utils.run_ppanggolin import run_ppanggolin_command

logger = logging.getLogger(__name__)


# Original ppanggolin commands run in this test used to be launched in CI main.yaml as:
"""
ppanggolin annotate --fasta genomes.fasta.list --output stepbystep --kingdom bacteria --cpu $NUM_CPUS
ppanggolin cluster -p stepbystep/pangenome.h5 --coverage 0.8 --identity 0.8 --cpu $NUM_CPUS
ppanggolin graph -p stepbystep/pangenome.h5 -r 10
ppanggolin partition --output stepbystep -f -p stepbystep/pangenome.h5 --cpu $NUM_CPUS -b 2.6 -ms 10 -fd -ck 500 -Kmm 3 12 -im 0.04 --draw_ICL

ppanggolin rgp -p stepbystep/pangenome.h5 --persistent_penalty 2 --variable_gain 1 --min_score 3 --dup_margin 0.05
ppanggolin spot -p stepbystep/pangenome.h5 --output stepbystep --spot_graph --overlapping_match 2 --set_size 3 --exact_match_size 1 -f
ppanggolin module -p stepbystep/pangenome.h5 --transitive 4 --size 3 --jaccard 0.86 --dup_margin 0.05

ppanggolin metrics -p stepbystep/pangenome.h5 --genome_fluidity --no_print_info --recompute_metrics --log metrics.log
ppanggolin info --pangenome stepbystep/pangenome.h5 > info_to_test/stepbystep_info.yaml
cat info_to_test/stepbystep_info.yaml
gzip -d stepbystep/gene_families.tsv.gz
echo "$(grep 'stepbystep/gene_families.tsv' expected_info_files/checksum.txt | cut -d' ' -f1)  stepbystep/gene_families.tsv" | shasum -a 256 -c - || { echo 'Checksum verification failed.' >&2; exit 1; }
shasum -a 256 stepbystep/gene_families.tsv >> info_to_test/checksum.txt
cd -

"""


@pytest.fixture(scope="session")
def stepbystep_pangenome(tmp_path_factory, num_cpus, request):
    """Run ppanggolin step-by-step workflow once and cache the result across pytest runs."""
    cache = request.config.cache
    cached_path = cache.get("ppanggolin/stepbystep_dir", None)

    if cached_path and Path(cached_path).exists():
        outdir = Path(cached_path)
        logger.warning("Reusing cached stepbystep workflow at %s", outdir)
    else:
        outdir = tmp_path_factory.mktemp("stepbystep_test") / "stepbystep"
        pangenome = outdir / "pangenome.h5"

        commands = [
            f"ppanggolin annotate --fasta testingDataset/genomes.fasta.list "
            f"--output {outdir} --kingdom bacteria --cpu {num_cpus}",
            f"ppanggolin cluster -p {pangenome} --coverage 0.8 --identity 0.8 --cpu {num_cpus}",
            f"ppanggolin graph -p {pangenome} -r 10",
            f"ppanggolin partition --output {outdir} -f -p {pangenome} --cpu {num_cpus} "
            f"-b 2.6 -ms 10 -fd -ck 500 -Kmm 3 12 -im 0.04 --draw_ICL",
            f"ppanggolin rgp -p {pangenome} --persistent_penalty 2 --variable_gain 1 "
            f"--min_score 3 --dup_margin 0.05",
            f"ppanggolin spot -p {pangenome} --output {outdir} --spot_graph "
            f"--overlapping_match 2 --set_size 3 --exact_match_size 1 -f",
            f"ppanggolin module -p {pangenome} --transitive 4 --size 3 "
            f"--jaccard 0.86 --dup_margin 0.05",
            f"ppanggolin metrics -p {pangenome} --genome_fluidity --no_print_info "
            f"--recompute_metrics --log metrics.log",
        ]

        for cmd in commands:
            logger.info("Running step-by-step command: %s", cmd)
            run_ppanggolin_command(cmd)

        cache.set("ppanggolin/stepbystep_dir", str(outdir))
        logger.info(
            f"Cached step-by-step workflow output at {outdir}",
        )

    return outdir / "pangenome.h5"


def test_gene_families(stepbystep_pangenome, update_golden):

    pan_dir = stepbystep_pangenome.parent
    cmds = [
        f"ppanggolin write_pangenome -p {stepbystep_pangenome} --output {pan_dir} -f --families_tsv --compress",
    ]
    for cmd in cmds:
        run_ppanggolin_command(cmd)

    file_path = pan_dir / "gene_families.tsv.gz"
    assert file_path.exists(), f"{file_path} does not exist"
    assert_or_update_hash(file_path, update_golden)


def test_spot_output(stepbystep_pangenome):

    pangenome_dir = stepbystep_pangenome.parent
    assert stepbystep_pangenome.exists()

    assert (pangenome_dir / "spotGraph.gexf").exists()
    assert (pangenome_dir / "spotGraph.gexf").stat().st_size > 0


def test_info_command(stepbystep_pangenome, tmp_path, update_golden):

    info_file = tmp_path / "stepbystep_info.yaml"

    check_pangenome_info(
        pangenome=stepbystep_pangenome, info_file=info_file, update_golden=update_golden
    )
