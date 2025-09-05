# tests/functional_tests/test_stepbystep.py
from pathlib import Path

from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_stepbystep import stepbystep_pangenome

from tests.functional_tests.test_fasta_all_wf import fasta_all_wf_pangenome
from tests.functional_tests.test_gbff_basic_wf import gbff_wf_pangenome


"""
ppanggolin write_genomes  -p stepbystep/pangenome.h5 --output stepbystep -f --fasta genomes.fasta.list --gff --proksee --table


head genomes.gbff.list | cut -f1  > genome_names.gbff.head.list

ppanggolin write_genomes  -p myannopang/pangenome.h5 --output flat_genomes_from_genome_files -f \
                            --anno genomes.gbff.list --gff --table --genomes  genome_names.gbff.head.list 

ppanggolin write_genomes  -p stepbystep/pangenome.h5 --output flat_genomes_from_cmdline_genomes --proksee \
                        --genomes GCF_006508185.1_ASM650818v1_genomic,GCF_002088315.1_ASM208831v1_genomic

head genomes.fasta.list | cut -f1  > genome_names.fasta.head.list 
# Default separator is a pipe but a pipe is found in a value of metadata db1. That is why we use another separator here. 
ppanggolin write_genomes -p mybasicpangenome/pangenome.h5 --output mybasicpangenome/genomes_outputs \
                        --genomes genome_names.fasta.head.list \
                            -f --gff --add_metadata --table --metadata_sep § --proksee

# Pipe separatore is found in metadata source db1. if we don't require this source then the writting with pipe is work fine. 
ppanggolin write_genomes -p mybasicpangenome/pangenome.h5 --output mybasicpangenome/genomes_outputs_with_metadata -f --gff --proksee --table --add_metadata  --metadata_sources db2 db3 db4 


"""


def _check_subdirs(outdir: Path, expected: dict):
    """Helper to check subdirectories and file counts."""
    for sub, count in expected.items():
        subdir = outdir / sub
        assert subdir.exists() and subdir.is_dir(), f"Missing subdir {subdir}"
        files = [f for f in subdir.iterdir() if f.is_file()]
        assert (
            len(files) == count
        ), f"Expected {count} files in {subdir}, found {len(files)}"


def test_write_genomes(stepbystep_pangenome, tmp_path):
    outdir = tmp_path / "write_genomes_outdir"
    cmd = (
        f"ppanggolin write_genomes -p {stepbystep_pangenome} "
        f"--output {outdir} -f --fasta testingDataset/genomes.fasta.list "
        f"--gff --proksee --table"
    )
    run_ppanggolin_command(cmd)
    _check_subdirs(outdir, {"gff": 51, "proksee": 51, "table": 51})


def test_write_genomes_from_anno_list(gbff_wf_pangenome, tmp_path):
    genome_list = tmp_path / "genome_names.gbff.head.list"
    # emulate: head genomes.gbff.list | cut -f1
    with (
        open("testingDataset/genomes.gbff.list") as f_in,
        open(genome_list, "w") as f_out,
    ):
        for i, line in enumerate(f_in):
            if i >= 10:  # limit to first 10 for test speed
                break
            f_out.write(line.split()[0] + "\n")

    outdir = tmp_path / "flat_genomes_from_genome_files"
    cmd = (
        f"ppanggolin write_genomes -p {gbff_wf_pangenome} --output {outdir} -f "
        f"--anno testingDataset/genomes.gbff.list --gff --table --genomes {genome_list}"
    )
    run_ppanggolin_command(cmd)
    _check_subdirs(outdir, {"gff": 10, "table": 10})


def test_write_genomes_from_cmdline(stepbystep_pangenome, tmp_path):
    outdir = tmp_path / "flat_genomes_from_cmdline_genomes"
    genomes = "GCF_006508185.1_ASM650818v1_genomic,GCF_002088315.1_ASM208831v1_genomic"
    cmd = (
        f"ppanggolin write_genomes -p {stepbystep_pangenome} --output {outdir} "
        f"--proksee --genomes {genomes}"
    )
    run_ppanggolin_command(cmd)
    _check_subdirs(outdir, {"proksee": 2})


def test_write_genomes_with_custom_metadata_sep(fasta_all_wf_pangenome, tmp_path):
    genome_list = tmp_path / "genome_names.fasta.head.list"
    with (
        open("testingDataset/genomes.fasta.list") as f_in,
        open(genome_list, "w") as f_out,
    ):
        for i, line in enumerate(f_in):
            if i >= 5:  # only first 5 genomes
                break
            f_out.write(line.split()[0] + "\n")

    outdir = tmp_path / "genomes_outputs"
    cmd = (
        f"ppanggolin write_genomes -p {fasta_all_wf_pangenome} --output {outdir} "
        f"--genomes {genome_list} -f --gff --add_metadata --table "
        f"--metadata_sep § --proksee"
    )
    run_ppanggolin_command(cmd)
    _check_subdirs(outdir, {"gff": 5, "table": 5, "proksee": 5})


def test_write_genomes_with_metadata_sources(fasta_all_wf_pangenome, tmp_path):
    outdir = tmp_path / "genomes_outputs_with_metadata"
    cmd = (
        f"ppanggolin write_genomes -p {fasta_all_wf_pangenome} --output {outdir} "
        f"-f --gff --proksee --table --add_metadata --metadata_sources db2 db3 db4"
    )
    run_ppanggolin_command(cmd)
    # we don’t know the exact number here, just ensure dirs exist and non-empty
    _check_subdirs(outdir, {"gff": 51, "proksee": 51, "table": 51})
