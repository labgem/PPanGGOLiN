import pytest
from pathlib import Path
from ppanggolin.annotate.annotate import extract_positions, read_anno_file


@pytest.mark.parametrize("input_string, expected_positions, expected_complement, expected_pseudogene", [
    ("join(190..7695,7695..12071)", [(190, 7695), (7695, 12071)], False, False),
    ("complement(join(4359800..4360707,4360707..4360962,1..100))", [(4359800, 4360707), (4360707, 4360962), (1, 100)],
     True, False),
    ("join(6835405..6835731,1..1218)", [(6835405, 6835731), (1, 1218)], False, False),
    ("join(1375484..1375555,1375557..1376579)", [(1375484, 1375555), (1375557, 1376579)], False, False),
    ("complement(6815492..6816265)", [(6815492, 6816265)], True, False),
    ("6811501..6812109", [(6811501, 6812109)], False, False),
    ("complement(6792573..>6795461)", [(6792573, 6795461)], True, True),
    ("join(1038313,1..1016)", [(1038313, 1038313), (1, 1016)], False, False),
    ("1038313", [(1038313, 1038313)], False, False)
])
def test_extract_positions(input_string, expected_positions, expected_complement, expected_pseudogene):
    positions, is_complement, is_pseudo = extract_positions(input_string)
    assert positions == expected_positions
    assert is_complement == expected_complement
    assert is_pseudo == expected_pseudogene


def test_extract_positions_with_wrong_positions_format():
    with pytest.raises(ValueError):
        extract_positions("join(1038313,1..1016")  # string misses a closing parenthesis


def test_extract_positions_with_wrong_positions_format2():
    with pytest.raises(ValueError):
        extract_positions("start..stop")  # start and stop are not integer
    with pytest.raises(ValueError):
        extract_positions("complement(join(start..6816265, 1..stop))")  # start and stop are not integer


@pytest.fixture
def genome_data():
    """
    Fixture providing common data for the tests.
    """
    # from the dir of the current python file, find a gbff file stored in testingDataset/GBFF repository
    genome_path = Path(__file__).resolve().parent / "../../testingDataset/GBFF/GCF_000026905.1_ASM2690v1_genomic.gbff.gz"
    circular_contigs = []
    genome_name = "GCF_000026905"
    return genome_name, genome_path, circular_contigs


def test_read_anno_file(genome_data):
    """
    Test reading annotation file without considering pseudogenes.
    """
    genome_name, genome_path, circular_contigs = genome_data
    use_pseudogene = False

    genome, has_sequence = read_anno_file(genome_name, genome_path, circular_contigs, use_pseudogene)

    assert has_sequence is True
    assert genome.name == genome_name
    assert genome.number_of_genes() == 905
    assert genome.number_of_contigs == 1


def test_read_anno_file_with_pseudo_enable(genome_data):
    """
    Test reading annotation file considering pseudogenes.
    """
    genome_name, genome_path, circular_contigs = genome_data
    use_pseudogene = True

    genome, has_sequence = read_anno_file(genome_name, genome_path, circular_contigs, use_pseudogene)

    assert has_sequence is True
    assert genome.name == genome_name
    assert genome.number_of_genes() == 928
    assert genome.number_of_contigs == 1


def test_with_joined_genes():
    genome_path = Path(__file__).resolve().parent / "../../testingDataset/GBFF/GCF_002776845.1_ASM277684v1_genomic.gbff.gz"
    genome_name = "GCF_002776845"
    circular_contigs = []
    use_pseudogene = True
    genome, _ = read_anno_file(genome_name, genome_path, circular_contigs, use_pseudogene)

    # this genome has 2 genes that are joined. 
    assert genome.number_of_genes() == 917
