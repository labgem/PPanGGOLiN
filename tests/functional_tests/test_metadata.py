# tests/functional_tests/test_stepbystep.py
from pathlib import Path
import pytest
import shutil
from tests.utils.run_ppanggolin import run_ppanggolin_command
from tests.functional_tests.test_fasta_all_wf import fasta_all_wf_pangenome

"""
ppanggolin metadata -p mybasicpangenome/pangenome.h5 -s db1 -m metadata/metadata_genes.tsv -a genes
ppanggolin metadata -p mybasicpangenome/pangenome.h5 -s db2 -m metadata/metadata_genomes.tsv -a genomes
ppanggolin metadata -p mybasicpangenome/pangenome.h5 -s db3 -m metadata/metadata_families.tsv -a families --omit
ppanggolin metadata -p mybasicpangenome/pangenome.h5 -s db4 -m metadata/metadata_rgps.tsv -a RGPs
ppanggolin metadata -p mybasicpangenome/pangenome.h5 -s db5 -m metadata/metadata_contigs.tsv  -a contigs
ppanggolin metadata -p mybasicpangenome/pangenome.h5 -s db6 -m metadata/metadata_modules.tsv  -a modules
ppanggolin write_metadata -p mybasicpangenome/pangenome.h5 -o metadata_flat_output


ppanggolin write_pangenome -p mybasicpangenome/pangenome.h5 --output mybasicpangenome -f --gexf --light_gexf --cpu $NUM_CPUS
ppanggolin rgp_cluster --pangenome mybasicpangenome/pangenome.h5 -o rgp_cluster_with_metadata --graph_formats graphml

cd -
"""


METADATA_DB_AND_FILES = [
    ("db1", "testingDataset/metadata/metadata_genes.tsv", "genes"),
    ("db2", "testingDataset/metadata/metadata_genomes.tsv", "genomes"),
    ("db3", "testingDataset/metadata/metadata_families.tsv", "families"),
    ("db4", "testingDataset/metadata/metadata_rgps.tsv", "RGPs"),
    ("db5", "testingDataset/metadata/metadata_contigs.tsv", "contigs"),
    ("db6", "testingDataset/metadata/metadata_modules.tsv", "modules"),
]

@pytest.fixture(scope="session")
def pangenome_with_all_metadata(fasta_all_wf_pangenome, tmp_path_factory):
    """
    Add all metadata to a fresh pangenome and return the modified HDF5 file.
    """
    tmp_path = tmp_path_factory.mktemp("metadata_test")
    fresh_outdir = tmp_path / "fresh_pangenome"
    fresh_outdir.mkdir(parents=True, exist_ok=True)

    # copy the pangenome HDF5
    pangenome = fresh_outdir / fasta_all_wf_pangenome.name
    shutil.copy2(fasta_all_wf_pangenome, pangenome)


    outdir = tmp_path / "metadata_added"
    outdir.mkdir(exist_ok=True)

    for db, metadata_file, attribute in METADATA_DB_AND_FILES:
        cmd = f"ppanggolin metadata -p {pangenome} -s {db} -m {metadata_file} -a {attribute}"
        if db == "db3":  # families require --omit
            cmd += " --omit"
        run_ppanggolin_command(cmd)

    return pangenome

def test_write_metadata(pangenome_with_all_metadata, tmp_path):
    """
    Test ppanggolin write_metadata command to produce flat metadata files.
    """
    outdir = tmp_path / "metadata_flat_output"

    cmd = f"ppanggolin write_metadata -p {pangenome_with_all_metadata} -o {outdir}"
    run_ppanggolin_command(cmd)

    # Check that some expected metadata files exist
    expected_files = [
        "contigs_metadata_from_db5.tsv",
        "families_metadata_from_db3.tsv",
        "genes_metadata_from_db1.tsv",
        "genomes_metadata_from_db2.tsv",
        "modules_metadata_from_db6.tsv",
        "RGPs_metadata_from_db4.tsv",
    ]

    for fname in expected_files:
        fpath = outdir / fname
        assert fpath.exists(), f"Expected metadata file {fname} not found"
        assert fpath.stat().st_size > 0, f"Metadata file {fname} is empty"
