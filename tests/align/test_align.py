import pytest
from typing import List
from random import choice, randint

from ppanggolin.align.alignOnPang import get_seq_ids


@pytest.fixture
def nucleotides():
    yield ['A', 'C', 'G', 'T', 'N']

@pytest.fixture
def aminoacids():
    amino_acid_alphabet = [
        'A',  # Alanine
        'R',  # Arginine
        'N',  # Asparagine
        'D',  # Aspartic acid
        'C',  # Cysteine
        'E',  # Glutamic acid
        'Q',  # Glutamine
        'G',  # Glycine
        'H',  # Histidine
        'I',  # Isoleucine
        'L',  # Leucine
        'K',  # Lysine
        'M',  # Methionine
        'F',  # Phenylalanine
        'P',  # Proline
        'S',  # Serine
        'T',  # Threonine
        'W',  # Tryptophan
        'Y',  # Tyrosine
        'V'  # Valine
    ]
    yield amino_acid_alphabet

@pytest.fixture
def number_of_sequences():
    yield randint(4, 10)
@pytest.fixture
def single_line_fasta(request, tmp_path_factory: pytest.TempPathFactory, number_of_sequences,
                      aminoacids: List[str], nucleotides: List[str]):
    if request.node.get_closest_marker('aminoacids'):
        alphabet = aminoacids
        fasta_path = tmp_path_factory.getbasetemp() / "single_line_nt.fasta"
    else:
        alphabet = nucleotides
        fasta_path = tmp_path_factory.getbasetemp() / "single_line_aa.fasta"
    with open(fasta_path, "w") as fasta:
        for i in range(number_of_sequences):
            fasta.write(f">Gene_{i}\n")
            fasta.write("".join([choice(alphabet) for _ in range(0, randint(30, 100))]))
            fasta.write("\n")
    yield fasta_path

@pytest.fixture
def multi_line_fasta(request, tmp_path_factory: pytest.TempPathFactory, number_of_sequences,
                     aminoacids: List[str], nucleotides: List[str]):
    if request.node.get_closest_marker('aminoacids'):
        alphabet = aminoacids
        fasta_path = tmp_path_factory.getbasetemp() / "single_line_nt.fasta"
    else:
        alphabet = nucleotides
        fasta_path = tmp_path_factory.getbasetemp() / "single_line_aa.fasta"
    with open(fasta_path, "w") as fasta:
        for i in range(number_of_sequences):
            fasta.write(f">Gene_{i}\n")
            for j in range(randint(4,10)):
                fasta.write("".join([choice(alphabet) for _ in range(60)]))
                fasta.write("\n")
    yield fasta_path

@pytest.mark.nucleotides
def test_get_seq_ids_single_line_nt(number_of_sequences, single_line_fasta):
    with open(single_line_fasta, "r") as fasta:
        seq_set, is_nucleotide, single_line_fasta = get_seq_ids(fasta)
        assert len(seq_set) == number_of_sequences
        assert seq_set == {f"Gene_{i}" for i in range(number_of_sequences)}
        assert is_nucleotide
        assert single_line_fasta

@pytest.mark.aminoacids
def test_get_seq_ids_single_line_aa(number_of_sequences, single_line_fasta):
    with open(single_line_fasta, "r") as fasta:
        seq_set, is_nucleotide, single_line_fasta = get_seq_ids(fasta)
        assert len(seq_set) == number_of_sequences
        assert seq_set == {f"Gene_{i}" for i in range(number_of_sequences)}
        assert not is_nucleotide
        assert single_line_fasta

@pytest.mark.nucleotides
def test_get_seq_ids_multi_line_nt(number_of_sequences, multi_line_fasta):
    with open(multi_line_fasta, "r") as fasta:
        seq_set, is_nucleotide, single_line_fasta = get_seq_ids(fasta)
        assert len(seq_set) == number_of_sequences
        assert seq_set == {f"Gene_{i}" for i in range(number_of_sequences)}
        assert is_nucleotide
        assert not single_line_fasta

@pytest.mark.aminoacids
def test_get_seq_ids_multi_line_aa(number_of_sequences, multi_line_fasta):
    with open(multi_line_fasta, "r") as fasta:
        seq_set, is_nucleotide, single_line_fasta = get_seq_ids(fasta)
        assert len(seq_set) == number_of_sequences
        assert seq_set == {f"Gene_{i}" for i in range(number_of_sequences)}
        assert not is_nucleotide
        assert not single_line_fasta



