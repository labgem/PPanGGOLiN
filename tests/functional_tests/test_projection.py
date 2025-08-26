# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import pytest

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_fasta_all_wf import fasta_all_wf_pangenome
from tests.functional_tests.test_gbff_basic_wf import gbff_wf_pangenome


"""
head genomes.gbff.list | sed 's/^/input_genome_/g' > genomes.gbff.head.list
ppanggolin projection --pangenome stepbystep/pangenome.h5  -o projection_from_list_of_gbff --anno genomes.gbff.head.list --gff --proksee --cpu $NUM_CPUS

head genomes.fasta.list | sed 's/^/input_genome_/g' > genomes.fasta.head.list
ppanggolin projection --pangenome myannopang/pangenome.h5  -o projection_from_list_of_fasta --fasta genomes.fasta.head.list --gff --proksee --cpu $NUM_CPUS

ppanggolin projection --pangenome mybasicpangenome/pangenome.h5  -o projection_from_single_fasta \
                        --genome_name chlam_A --fasta FASTA/GCF_002776845.1_ASM277684v1_genomic.fna.gz \
                        --spot_graph --graph_formats graphml --fast --keep_tmp -f --add_sequences --gff --proksee --table --add_metadata --cpu $NUM_CPUS

ppanggolin projection --pangenome mybasicpangenome/pangenome.h5  -o projection_from_gff_prodigal \
                        --genome_name chlam_annotated_with_prodigal --anno GBFF/GCF_003788785.1_ct114V1_genomic_prodigal_annotation.gff.gz \
                        --gff  --table --cpu $NUM_CPUS

# projection of a plasmid with chevron that have been added manually to test chevron handeling in GFF
ppanggolin projection --pangenome myannopang/pangenome.h5 --anno GBFF/plasmid_NZ_CP007132_with_manually_added_chevrons.gff.gz --cpu $NUM_CPUS -o projection_plasmid_with_chevron

# projection with GFF with no sequence and fasta sequence
ppanggolin projection -p myannopang/pangenome.h5 --anno GBFF/plasmid_GCF_000093005.1_ASM9300v1.gff.gz --fasta GBFF/plasmid_GCF_000093005.1_ASM9300v1.fna.gz

# projection with GFF with no sequence and fasta sequence specified in a TSV file with other GFF (but with sequences)
head -n 3 genomes.gbff.head.list > genomes.gbff.h3_and_GFFplasmidNoSeq.list

echo GFF_plasmid_No_seq$'\t'GBFF/plasmid_GCF_000093005.1_ASM9300v1.gff.gz >> genomes.gbff.h3_and_GFFplasmidNoSeq.list
echo GFF_plasmid_No_seq$'\t'GBFF/plasmid_GCF_000093005.1_ASM9300v1.fna.gz >> genomes.fna.GFFplasmidNoSeq.list
ppanggolin projection -p myannopang/pangenome.h5 --anno genomes.gbff.h3_and_GFFplasmidNoSeq.list --fasta  genomes.fna.GFFplasmidNoSeq.list


# test projection with GBK from pyrodigal. gbk from pyrodigal was producing an error on the annotation step see issue #339 
ppanggolin projection -p myannopang/pangenome.h5 --anno GBFF/plasmid_NC_017433.1_pyrodigal.gbk.gz 
"""


@pytest.fixture
def make_head_list(tmp_path):
    def _make_head_list(
        original_file: str,
        prefix: str = "input_genome_",
        n: int = 10,
        data_dir="testingDataset",
    ) -> Path:
        """Create a head list file with a prefix (like `head file | sed 's/^/prefix/g'`)."""
        head_file = tmp_path / (Path(original_file).name + ".head.list")
        with open(original_file) as fin, open(head_file, "w") as fout:
            for i, line in enumerate(fin):
                if i >= n or not line:
                    break

                genome_name, path, *rest = line.strip().split("\t")
                line = f"{prefix}{genome_name}\t{data_dir}/{path}\n"
                fout.write(line)
        return head_file

    return _make_head_list


def test_projection_from_list_of_gbff(
    fasta_all_wf_pangenome, make_head_list, tmp_path, num_cpus
):
    head_file = make_head_list("testingDataset/genomes.gbff.list")
    outdir = tmp_path / "projection_from_list_of_gbff"

    cmd = (
        f"ppanggolin projection --pangenome {fasta_all_wf_pangenome} "
        f"-o {outdir} --anno {head_file} --gff --proksee --cpu {num_cpus}"
    )
    run_ppanggolin_command(cmd)

    assert outdir.exists()


def test_projection_from_list_of_fasta(
    gbff_wf_pangenome, make_head_list, tmp_path, num_cpus
):
    head_file = make_head_list("testingDataset/genomes.fasta.list")
    outdir = tmp_path / "projection_from_list_of_fasta"

    cmd = (
        f"ppanggolin projection --pangenome {gbff_wf_pangenome} "
        f"-o {outdir} --fasta {head_file} --gff --proksee --cpu {num_cpus}"
    )
    run_ppanggolin_command(cmd)

    assert outdir.exists()


def test_projection_from_single_fasta(fasta_all_wf_pangenome, tmp_path, num_cpus):
    outdir = tmp_path / "projection_from_single_fasta"
    cmd = (
        f"ppanggolin projection --pangenome {fasta_all_wf_pangenome} -o {outdir} "
        "--genome_name chlam_A --fasta testingDataset/FASTA/GCF_002776845.1_ASM277684v1_genomic.fna.gz "
        "--spot_graph --graph_formats graphml --fast --keep_tmp -f "
        "--add_sequences --gff --proksee --table --add_metadata --cpu {num_cpus}"
    ).format(num_cpus=num_cpus)

    run_ppanggolin_command(cmd)

    assert outdir.exists()


def test_projection_from_gff_prodigal(fasta_all_wf_pangenome, tmp_path, num_cpus):
    outdir = tmp_path / "projection_from_gff_prodigal"
    cmd = (
        f"ppanggolin projection --pangenome {fasta_all_wf_pangenome} -o {outdir} "
        "--genome_name chlam_annotated_with_prodigal "
        "--anno testingDataset/GBFF/GCF_003788785.1_ct114V1_genomic_prodigal_annotation.gff.gz "
        "--gff --table --cpu {num_cpus}"
    ).format(num_cpus=num_cpus)

    run_ppanggolin_command(cmd)
    assert outdir.exists()


def test_projection_plasmid_with_chevron(gbff_wf_pangenome, tmp_path, num_cpus):
    outdir = tmp_path / "projection_plasmid_with_chevron"
    cmd = (
        f"ppanggolin projection --pangenome {gbff_wf_pangenome} "
        "--anno testingDataset/GBFF/plasmid_NZ_CP007132_with_manually_added_chevrons.gff.gz "
        f"--cpu {num_cpus} -o {outdir}"
    )
    run_ppanggolin_command(cmd)


def test_projection_gff_plus_fasta(gbff_wf_pangenome, tmp_path):
    outdir = tmp_path / "projection_gff_plus_fasta"
    cmd = (
        f"ppanggolin projection -p {gbff_wf_pangenome} "
        "--anno testingDataset/GBFF/plasmid_GCF_000093005.1_ASM9300v1.gff.gz "
        "--fasta testingDataset/GBFF/plasmid_GCF_000093005.1_ASM9300v1.fna.gz "
        f"-o {outdir}"
    )
    run_ppanggolin_command(cmd)


def test_projection_from_list_of_gbff(
    fasta_all_wf_pangenome, make_head_list, tmp_path, num_cpus
):

    gbff_head_file = make_head_list("testingDataset/genomes.gbff.list", n=3)
    fna_file = tmp_path / "genomes.fna.GFFplasmidNoSeq.list"
    plasmide_file_basename = "testingDataset/GBFF/plasmid_GCF_000093005.1_ASM9300v1"
    plasmide_name = "GFF_plasmid_No_seq"

    with open(gbff_head_file, "a") as fout:
        fout.write(f"{plasmide_name}\t{plasmide_file_basename}.gff.gz\n")

    with open(fna_file, "w") as fout:
        fout.write(f"{plasmide_name}\t{plasmide_file_basename}.fna.gz\n")

    outdir = tmp_path / "projection_from_list_of_gbff"

    outdir = tmp_path / "projection_gff_tsv_plus_fasta_tsv"
    cmd = (
        f"ppanggolin projection -p {fasta_all_wf_pangenome} "
        f"--anno {gbff_head_file} "
        f"--fasta {fna_file} "
        f"-o {outdir} --cpu {num_cpus}"
    )
    run_ppanggolin_command(cmd)


def test_projection_pyrodigal_gbk(gbff_wf_pangenome, tmp_path):
    outdir = tmp_path / "projection_pyrodigal_gbk"
    cmd = (
        f"ppanggolin projection -p {gbff_wf_pangenome} "
        "--anno testingDataset/GBFF/plasmid_NC_017433.1_pyrodigal.gbk.gz "
        f"-o {outdir}"
    )
    run_ppanggolin_command(cmd)
