import pytest
from typing import List
from random import choice, randint

from ppanggolin.align.alignOnPang import get_seq_ids


@pytest.fixture
def single_line_fasta_nt() -> List:

    return [
        ">Gene_1 seq_description\n",
        "ATGCGTTGTCGTTG\n",
        ">Gene_2\n",
        "TGTGACCTGCT\n",
    ]


@pytest.fixture
def single_line_fasta_aa() -> List:

    return [
        ">Gene_1 seq_description\n",
        "YWTPRPFFYAAEYNN\n",
        ">Gene_2\n",
        "YWTPRPSYWTPAAEYNN\n",
    ]


@pytest.fixture
def multi_line_fasta_nt() -> List:

    return [
        ">Gene_1 seq_description\n",
        "ATGCGT\n",
        "TGTCGTTG\n",
        ">Gene_2\n",
        "TGTGACCTGCT\n",
    ]


@pytest.fixture
def multi_line_fasta_aa() -> List:

    return [
        ">Gene_1 seq_description\n",
        "AAEYNN\n",
        "YWTPRPFFY\n",
        ">Gene_2\n",
        "YWTPRPS\n",
        "YWTPAAEYNN\n",
    ]


def test_get_seq_ids_single_line_nt(single_line_fasta_nt):

    seq_set, is_nucleotide, single_line_fasta = get_seq_ids(single_line_fasta_nt)
    assert len(seq_set) == 2
    assert seq_set == {"Gene_1", "Gene_2"}
    assert is_nucleotide
    assert single_line_fasta


def test_get_seq_ids_single_line_aa(single_line_fasta_aa):
    seq_set, is_nucleotide, single_line_fasta = get_seq_ids(single_line_fasta_aa)
    assert len(seq_set) == 2
    assert seq_set == {"Gene_1", "Gene_2"}
    assert not is_nucleotide
    assert single_line_fasta


def test_get_seq_ids_multi_line_nt(multi_line_fasta_nt):

    seq_set, is_nucleotide, single_line_fasta = get_seq_ids(multi_line_fasta_nt)

    assert len(seq_set) == 2
    assert seq_set == {"Gene_1", "Gene_2"}
    assert is_nucleotide
    assert not single_line_fasta


def test_get_seq_ids_multi_line_aa(multi_line_fasta_aa):

    seq_set, is_nucleotide, single_line_fasta = get_seq_ids(multi_line_fasta_aa)
    assert len(seq_set) == 2
    assert seq_set == {"Gene_1", "Gene_2"}
    assert not is_nucleotide
    assert not single_line_fasta
