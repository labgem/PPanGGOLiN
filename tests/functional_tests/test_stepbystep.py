# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import pytest

from tests.utils.checksum import assert_or_update_hash
from tests.utils.file_compare import check_pangenome_info
from tests.utils.run_ppanggolin import run_ppanggolin_command


# Original ppanggolin commands run in this test used to be launched in CI main.yaml as:
"""
- name: Step by Step workflow with most options calls
      shell: bash -l {0}
      run: |
        cd testingDataset
        ppanggolin annotate --fasta genomes.fasta.list --output stepbystep --kingdom bacteria --cpu $NUM_CPUS
        ppanggolin cluster -p stepbystep/pangenome.h5 --coverage 0.8 --identity 0.8 --cpu $NUM_CPUS
        ppanggolin graph -p stepbystep/pangenome.h5 -r 10
        ppanggolin partition --output stepbystep -f -p stepbystep/pangenome.h5 --cpu $NUM_CPUS -b 2.6 -ms 10 -fd -ck 500 -Kmm 3 12 -im 0.04 --draw_ICL
        ppanggolin rarefaction --output stepbystep -f -p stepbystep/pangenome.h5 --depth 5 --min 1 --max 50 -ms 10 -fd -ck 30 -K 3 --soft_core 0.9 -se $RANDOM

        # Check reestimate K with depth 1 and min 40 to avoid long computation time
        ppanggolin rarefaction --output stepbystep -f -p stepbystep/pangenome.h5 --reestimate_K --depth 1 --min 40

        ppanggolin draw -p stepbystep/pangenome.h5 --tile_plot --nocloud --soft_core 0.92 --ucurve --output stepbystep -f
        ppanggolin rgp -p stepbystep/pangenome.h5 --persistent_penalty 2 --variable_gain 1 --min_score 3 --dup_margin 0.05
        ppanggolin spot -p stepbystep/pangenome.h5 --output stepbystep --spot_graph --overlapping_match 2 --set_size 3 --exact_match_size 1 -f
        ppanggolin draw -p stepbystep/pangenome.h5 --draw_spots -o stepbystep -f
        ppanggolin module -p stepbystep/pangenome.h5 --transitive 4 --size 3 --jaccard 0.86 --dup_margin 0.05
        ppanggolin write_pangenome -p stepbystep/pangenome.h5 --output stepbystep -f --soft_core 0.9 --dup_margin 0.06  --gexf --light_gexf --csv --Rtab --stats --partitions --compress --json --spots --regions --borders --families_tsv --cpu 1 
        ppanggolin write_genomes  -p stepbystep/pangenome.h5 --output stepbystep -f --fasta genomes.fasta.list --gff --proksee --table
        ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --prot_families all --gene_families shell --regions all --fasta genomes.fasta.list
        ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --prot_families rgp --gene_families rgp --compress 
        ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --prot_families softcore --gene_families softcore 
        ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --prot_families module_0
        ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --genes core --proteins cloud
        ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --gene_families module_0 --genes module_0 --compress
        ppanggolin fasta -p stepbystep/pangenome.h5 --output stepbystep -f --proteins cloud --cpu $NUM_CPUS --keep_tmp --compress

        ppanggolin draw -p stepbystep/pangenome.h5 --draw_spots --spots all -o stepbystep -f
        ppanggolin metrics -p stepbystep/pangenome.h5 --genome_fluidity --no_print_info --recompute_metrics --log metrics.log
        ppanggolin info --pangenome stepbystep/pangenome.h5 > info_to_test/stepbystep_info.yaml
        cat info_to_test/stepbystep_info.yaml
        gzip -d stepbystep/gene_families.tsv.gz
        echo "$(grep 'stepbystep/gene_families.tsv' expected_info_files/checksum.txt | cut -d' ' -f1)  stepbystep/gene_families.tsv" | shasum -a 256 -c - || { echo 'Checksum verification failed.' >&2; exit 1; }
        shasum -a 256 stepbystep/gene_families.tsv >> info_to_test/checksum.txt
        cd -

"""

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
    run_ppanggolin_command(cmd)
    return outdir


@pytest.fixture(scope="module")
def rgp_and_spot(partition):
    pangenome = partition / "pangenome.h5"
    outdir = partition
    cmd = f"ppanggolin rgp -p {pangenome} --persistent_penalty 2 --variable_gain 1 --min_score 3 --dup_margin 0.05"
    run_ppanggolin_command(cmd)
    cmd = f"ppanggolin spot -p {pangenome} --output {outdir} --spot_graph --overlapping_match 2 --set_size 3 --exact_match_size 1 -f"
    run_ppanggolin_command(cmd)
    return outdir


@pytest.fixture(scope="module")
def module(rgp_and_spot):
    pangenome = rgp_and_spot / "pangenome.h5"
    outdir = rgp_and_spot
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


def test_rgp_and_rgp(rgp_and_spot):
    pangenome = rgp_and_spot / "pangenome.h5"
    assert pangenome.exists()

    assert (rgp_and_spot / "spotGraph.gexf").exists()
    assert (rgp_and_spot / "spotGraph.gexf").stat().st_size > 0


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


OUTDIR_COMMANDS_AND_FILES = [
    (
        "rarefaction1",
        "ppanggolin rarefaction --output {outdir} -f -p {pangenome_dir}/pangenome.h5 "
        "--depth 2 --min 1 --max 50 -ms 10 -fd -ck 30 -K 3 --soft_core 0.9 -se 42",
        [
            "rarefaction.csv",
            "rarefaction_curve.html",
            "rarefaction_parameters.csv",
        ],
    ),
    (
        "rarefaction2",
        "ppanggolin rarefaction --output {outdir} -f -p {pangenome_dir}/pangenome.h5 "
        "--reestimate_K --depth 1 --min 40",
        [
            "rarefaction.csv",
            "rarefaction_curve.html",
            "rarefaction_parameters.csv",
        ],
    ),
    (
        "draw_tile",
        "ppanggolin draw -p {pangenome_dir}/pangenome.h5 --tile_plot --nocloud --soft_core 0.92 "
        "--ucurve --output {outdir}",
        [
            "tile_plot.html",
            "Ushaped_plot.html",
        ],
    ),
    (
        "draw_spots",
        "ppanggolin draw -p {pangenome_dir}/pangenome.h5 --draw_spots -o {outdir} -f",
        [
            "spot_0.gexf",
            "spot_0.html",
            "spot_0_identical_rgps.tsv",
            "spot_1.gexf",
            "spot_1.html",
            "spot_1_identical_rgps.tsv",
        ],
    ),
    (
        "write_pangenome_outdir",
        "ppanggolin write_pangenome -p {pangenome_dir}/pangenome.h5 --output {outdir} -f "
        "--soft_core 0.9 --dup_margin 0.06 --gexf --light_gexf --csv --Rtab --stats "
        "--partitions --compress --json --spots --regions --borders --families_tsv --cpu 1",
        [
            "border_protein_genes.fasta.gz",
            "genomes_statistics.tsv.gz",
            "pangenomeGraph.gexf.gz",
            "gene_families.tsv.gz",
            "matrix.csv.gz",
            "pangenomeGraph.json.gz",
            "regions_of_genomic_plasticity.tsv.gz",
            "summarize_spots.tsv.gz",
            "gene_presence_absence.Rtab.gz",
            "mean_persistent_duplication.tsv.gz",
            "pangenomeGraph_light.gexf.gz",
            "spot_borders.tsv.gz",
            "spots.tsv.gz",
            "partitions/cloud.txt",
            "partitions/exact_accessory.txt",
            "partitions/exact_core.txt",
            "partitions/persistent.txt",
            "partitions/S1.txt",
            "partitions/shell.txt",
            "partitions/soft_accessory.txt",
            "partitions/soft_core.txt",
        ],
    ),
    (
        "fasta_all_shell",
        "ppanggolin fasta -p {pangenome_dir}/pangenome.h5 --output {outdir} -f "
        "--prot_families all --gene_families shell --regions all "
        "--fasta testingDataset/genomes.fasta.list",
        [
            "all_protein_families.faa",
            "all_rgp_genomic_sequences.fasta",
            "shell_nucleotide_families.fasta",
        ],
    ),
    (
        "fasta_rgp",
        "ppanggolin fasta -p {pangenome_dir}/pangenome.h5 --output {outdir} -f "
        "--prot_families rgp --gene_families rgp --compress",
        [
            "rgp_nucleotide_families.fasta.gz",
            "rgp_protein_families.faa.gz",
        ],
    ),
    (
        "fasta_softcore",
        "ppanggolin fasta -p {pangenome_dir}/pangenome.h5 --output {outdir} -f "
        "--prot_families softcore --gene_families softcore",
        [
            "softcore_nucleotide_families.fasta",
            "softcore_protein_families.faa",
        ],
    ),
    (
        "fasta_module0",
        "ppanggolin fasta -p {pangenome_dir}/pangenome.h5 --output {outdir} -f --prot_families module_0",
        [
            "module_0_protein_families.faa",
        ],
    ),
    (
        "fasta_core_cloud",
        "ppanggolin fasta -p {pangenome_dir}/pangenome.h5 --output {outdir} -f --genes core --proteins cloud",
        [
            "cloud_protein_genes.faa",
            "core_genes.fna",
        ],
    ),
    (
        "fasta_module0_compress",
        "ppanggolin fasta -p {pangenome_dir}/pangenome.h5 --output {outdir} -f "
        "--gene_families module_0 --genes module_0 --compress",
        [
            "module_0_genes.fna.gz",
            "module_0_nucleotide_families.fasta.gz",
        ],
    ),
    (
        "fasta_cloud_compress",
        "ppanggolin fasta -p {pangenome_dir}/pangenome.h5 --output {outdir} -f "
        "--proteins cloud --cpu 2 --keep_tmp --compress",
        [
            "cloud_protein_genes.faa.gz",
        ],
    ),
    (
        "draw_spots_all",
        "ppanggolin draw -p {pangenome_dir}/pangenome.h5 --draw_spots --spots all -o {outdir} -f",
        [
            "spot_0.gexf",
            "spot_0.html",
            "spot_0_identical_rgps.tsv",
            "spot_1.gexf",
            "spot_1.html",
            "spot_1_identical_rgps.tsv",
        ],
    ),
]


@pytest.mark.parametrize(
    "test_outdir, cmd_template, expected_files",
    OUTDIR_COMMANDS_AND_FILES,
    ids=[case[0] for case in OUTDIR_COMMANDS_AND_FILES],
)
def test_stepbystep_outputs(
    module, tmp_path, test_outdir, cmd_template, expected_files
):
    """Run each command on the prepared pangenome and check expected files exist."""
    outdir = tmp_path / test_outdir
    cmd = cmd_template.format(pangenome_dir=module, outdir=outdir)

    run_ppanggolin_command(cmd)

    for fname in expected_files:
        fpath = Path(outdir) / fname
        assert fpath.exists(), f"Expected file {fname} not found after `{cmd}`"
        assert fpath.stat().st_size > 0, f"File {fname} is empty after `{cmd}`"


def test_write_genomes(module):
    outdir = module / "write_genomes_outdir"
    pangenome = module / "pangenome.h5"
    cmd = (
        f"ppanggolin write_genomes "
        f"-p {pangenome} "
        f"--output {outdir} "
        f"-f --fasta testingDataset/genomes.fasta.list "
        f"--gff --proksee --table"
    )

    run_ppanggolin_command(cmd)

    # Check the expected subdirectories exist
    expected_dirs = ["gff", "proksee", "table"]
    for sub in expected_dirs:
        subdir = outdir / sub
        assert subdir.exists() and subdir.is_dir(), f"Missing subdir {subdir}"

        # Check number of files inside
        files = [f for f in subdir.iterdir() if f.is_file()]
        assert len(files) == 51, f"Expected 51 files in {subdir}, found {len(files)}"
