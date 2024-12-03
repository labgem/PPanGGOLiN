import pytest
from pathlib import Path

from ppanggolin.genome import Contig
from ppanggolin.annotate.annotate import (
    extract_positions,
    read_anno_file,
    parse_contig_header_lines,
    parse_gbff_by_contig,
    parse_feature_lines,
    parse_dna_seq_lines,
    read_org_gbff,
    combine_contigs_metadata,
    fix_partial_gene_coordinates,
    shift_start_coordinates,
    shift_end_coordinates,
)

from ppanggolin.annotate.synta import check_sequence_tuple, parse_fasta


@pytest.mark.parametrize(
    "input_string, expected_positions, expected_complement, expected_partialgene_start, expected_partialgene_end",
    [
        (
            "join(190..7695,7695..12071)",
            [(190, 7695), (7695, 12071)],
            False,
            False,
            False,
        ),
        (
            "order(190..7695,7995..12071)",
            [(190, 7695), (7995, 12071)],
            False,
            False,
            False,
        ),
        (
            "complement(join(4359800..4360707,4360707..4360962,1..100))",
            [(4359800, 4360707), (4360707, 4360962), (1, 100)],
            True,
            False,
            False,
        ),
        (
            "complement(order(4359800..4360707,4360707..4360962,1..100))",
            [(4359800, 4360707), (4360707, 4360962), (1, 100)],
            True,
            False,
            False,
        ),
        (
            "join(6835405..6835731,1..1218)",
            [(6835405, 6835731), (1, 1218)],
            False,
            False,
            False,
        ),
        (
            "join(1375484..1375555,1375557..1376579)",
            [(1375484, 1375555), (1375557, 1376579)],
            False,
            False,
            False,
        ),
        ("complement(6815492..6816265)", [(6815492, 6816265)], True, False, False),
        ("6811501..6812109", [(6811501, 6812109)], False, False, False),
        ("complement(6792573..>6795461)", [(6792573, 6795461)], True, False, True),
        ("complement(<6792573..6795461)", [(6792573, 6795461)], True, True, False),
        ("complement(<6792573..>6795461)", [(6792573, 6795461)], True, True, True),
        ("join(1038313,1..1016)", [(1038313, 1038313), (1, 1016)], False, False, False),
        ("1038313", [(1038313, 1038313)], False, False, False),
    ],
)
def test_extract_positions(
    input_string,
    expected_positions,
    expected_complement,
    expected_partialgene_start,
    expected_partialgene_end,
):
    positions, is_complement, has_partial_start, has_partial_end = extract_positions(
        input_string
    )
    assert positions == expected_positions
    assert is_complement == expected_complement
    assert has_partial_start == expected_partialgene_start
    assert has_partial_end == expected_partialgene_end


def test_extract_positions_with_wrong_positions_format():
    with pytest.raises(ValueError):
        extract_positions("join(1038313,1..1016")  # string misses a closing parenthesis


def test_extract_positions_with_strange_chevrons():
    with pytest.raises(ValueError):
        extract_positions(
            "complement(join(4359800..>4360707,1..100))"
        )  # chevron in inner position
    with pytest.raises(ValueError):
        extract_positions(
            "complement(join(4359800..4360707,<1..100))"
        )  # chevron in inner position

    with pytest.raises(ValueError):
        extract_positions(
            "complement(join(4359800..4360707,1..<100))"
        )  # start chevron in ending position


def test_extract_positions_with_wrong_positions_format2():
    with pytest.raises(ValueError):
        extract_positions("start..stop")  # start and stop are not integer
    with pytest.raises(ValueError):
        extract_positions(
            "complement(join(start..6816265, 1..stop))"
        )  # start and stop are not integer
    with pytest.raises(ValueError):
        extract_positions("start..stop")  # start and stop are not integer
    with pytest.raises(ValueError):
        extract_positions(
            "complement(join(start..6816265, 1..stop))"
        )  # start and stop are not integer


@pytest.fixture
def genome_data():
    """
    Fixture providing common data for the tests.
    """
    script_path = Path(__file__).resolve()
    ppanggolin_main_dir = script_path.parent.parent.parent

    genome_path = (
        ppanggolin_main_dir
        / "testingDataset/GBFF/GCF_000026905.1_ASM2690v1_genomic.gbff.gz"
    )
    circular_contigs = []
    genome_name = "GCF_000026905"
    return genome_name, genome_path, circular_contigs


@pytest.fixture
def genome_data_with_joined_genes():
    """
    Fixture providing gbff file with joined genes
    """
    script_path = Path(__file__).resolve()
    ppanggolin_main_dir = script_path.parent.parent.parent

    genome_path = (
        ppanggolin_main_dir
        / "testingDataset/GBFF/GCF_002776845.1_ASM277684v1_genomic.gbff.gz"
    )
    circular_contigs = []
    genome_name = "GCF_002776845"
    return genome_name, genome_path, circular_contigs


def test_read_anno_file(genome_data):
    """
    Test reading annotation file without considering pseudogenes.
    """
    genome_name, genome_path, circular_contigs = genome_data
    use_pseudogene = False

    genome, has_sequence = read_anno_file(
        genome_name, genome_path, circular_contigs, use_pseudogene
    )

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

    genome, has_sequence = read_anno_file(
        genome_name, genome_path, circular_contigs, use_pseudogene
    )

    assert has_sequence is True
    assert genome.name == genome_name
    assert genome.number_of_genes() == 928
    assert genome.number_of_contigs == 1


def test_with_joined_genes(genome_data_with_joined_genes):
    genome_name, genome_path, circular_contigs = genome_data_with_joined_genes
    use_pseudogene = True
    genome, _ = read_anno_file(
        genome_name, genome_path, circular_contigs, use_pseudogene
    )

    # this genome has 2 genes that are joined.
    assert genome.number_of_genes() == 917


def test_read_org_gbff(genome_data_with_joined_genes):
    genome_name, genome_path, circular_contigs = genome_data_with_joined_genes
    genome, _ = read_org_gbff(
        genome_name, genome_path, circular_contigs, use_pseudogenes=True
    )

    # this genome has 2 genes that are joined.
    assert genome.number_of_genes() == 917


def test_gbff_header_parser():
    header_lines = [
        "LOCUS       NC_022109            1041595 bp    DNA     circular CON 24-MAR-2017",
        "DEFINITION  Chlamydia trachomatis strain D/14-96 genome.",
        "VERSION     NC_022109.1",
        "SOURCE      Chlamydia trachomatis",
        "  ORGANISM  Chlamydia trachomatis",
        "            Bacteria; Chlamydiae; Chlamydiales; Chlamydiaceae;",
        "            Chlamydia/Chlamydophila group; Chlamydia.",
    ]

    parsed_header = parse_contig_header_lines(header_lines)

    assert parsed_header == {
        "LOCUS": "NC_022109            1041595 bp    DNA     circular CON 24-MAR-2017",
        "DEFINITION": "Chlamydia trachomatis strain D/14-96 genome.",
        "VERSION": "NC_022109.1",
        "SOURCE": "Chlamydia trachomatis",
        "ORGANISM": "Chlamydia trachomatis\nBacteria; Chlamydiae; Chlamydiales; Chlamydiaceae;\nChlamydia/Chlamydophila group; Chlamydia.",
    }


# Define test data
@pytest.fixture
def sample_gbff_path(tmp_path):
    gbff_content = """LOCUS       NC_022109            1041595 bp    DNA     circular CON 24-MAR-2017
DEFINITION  Chlamydia trachomatis strain D/14-96 genome.
ACCESSION   NC_022109
VERSION     NC_022109.1
KEYWORDS    RefSeq.
SOURCE      Chlamydia trachomatis
FEATURES             Location/Qualifiers
     source          1..1041595
                     /organism="Chlamydia trachomatis"
                     /mol_type="genomic DNA"
ORIGIN
        1 aaaccgggtt
       11 ccaaatttgg
//
LOCUS       NC_022110            2041595 bp    DNA     circular CON 24-MAR-2017
DEFINITION  Another genome.
KEYWORDS    RefSeq.
SOURCE      Chlamydia trachomatis
  ORGANISM  Chlamydia trachomatis
            Bacteria; Chlamydiae; Chlamydiales; Chlamydiaceae;
            Chlamydia/Chlamydophila group; Chlamydia.
FEATURES             Location/Qualifiers
     source          1..2041595
                     /organism="Chlamydia trachomatis"
                     /mol_type="genomic DNA"
ORIGIN
        1 aaaccgggtt
       11 ccaaatttgg
       21 ggcccctttt
//
"""

    gbff_file_path = tmp_path / "sample.gbff"
    with open(gbff_file_path, "w") as f:
        f.write(gbff_content)
    return gbff_file_path


# Define test cases
def test_parse_gbff_by_contig(sample_gbff_path):
    contigs = list(parse_gbff_by_contig(sample_gbff_path))

    assert len(contigs) == 2

    # Check first contig
    header_1, feature_1, sequence_1 = contigs[0]
    assert len(header_1) == 6
    assert len(list(feature_1)) == 1
    assert sequence_1 == "AAACCGGGTTCCAAATTTGG"

    # Check second contig
    header_2, feature_2, sequence_2 = contigs[1]
    assert len(header_2) == 5
    assert list(feature_2) == [
        {
            "feature_type": "source",
            "location": "1..2041595",
            "organism": "Chlamydia trachomatis",
            "mol_type": "genomic DNA",
        }
    ]

    assert sequence_2 == "AAACCGGGTTCCAAATTTGGGGCCCCTTTT"


# Define test cases
@pytest.mark.parametrize(
    "input_lines, expected_output",
    [
        (
            [
                "     gene            123..456",
                "                     /locus_tag=ABC123",
                "                     /note=Some note",
                "     CDS             789..1011",
                "                     /protein_id=DEF456",
                "                     /translation=ATGCTAGCATCG",
            ],
            [
                {
                    "feature_type": "gene",
                    "location": "123..456",
                    "locus_tag": "ABC123",
                    "note": "Some note",
                },
                {
                    "feature_type": "CDS",
                    "location": "789..1011",
                    "protein_id": "DEF456",
                    "translation": "ATGCTAGCATCG",
                },
            ],
        ),
        (
            [
                "     gene            123..456",
                "                     /locus_tag=ABC123",
                "                     /note=Some note",
                "     CDS             789..1011",
                '                     /protein_id="DEF456"',
                '                     /translation="ATGCTAGCATCG"',
                "     gene            789..1011",
                "                     /locus_tag=DEF789",
                "                     /note=Another note",
            ],
            [
                {
                    "feature_type": "gene",
                    "location": "123..456",
                    "locus_tag": "ABC123",
                    "note": "Some note",
                },
                {
                    "feature_type": "CDS",
                    "location": "789..1011",
                    "protein_id": "DEF456",
                    "translation": "ATGCTAGCATCG",
                },
                {
                    "feature_type": "gene",
                    "location": "789..1011",
                    "locus_tag": "DEF789",
                    "note": "Another note",
                },
            ],
        ),
        # Add more test cases as needed
    ],
)
def test_parse_feature_lines(input_lines, expected_output):
    assert list(parse_feature_lines(input_lines)) == expected_output


def test_parse_dna_seq_lines():
    lines = ["        1 aaacc gggtt", "       11 ccaaa tttgg", "       21 ggccc ctttt"]

    assert parse_dna_seq_lines(lines) == "AAACCGGGTTCCAAATTTGGGGCCCCTTTT"


def test_combine_contigs_metadata():
    contig1, contig2, contig3 = (
        Contig(1, "contig1"),
        Contig(2, "contig2"),
        Contig(3, "contig3"),
    )
    contig_to_metadata = {
        contig1: {"sp": "spA", "strain": "123", "contig_feat": "ABC"},
        contig2: {"sp": "spA", "strain": "123", "contig_feat": "XYZ"},
        contig3: {"sp": "spA", "strain": "123"},
    }

    genome_metadata, contig_metadata = combine_contigs_metadata(contig_to_metadata)

    assert genome_metadata == {"sp": "spA", "strain": "123"}
    assert contig_metadata == {
        contig1: {"contig_feat": "ABC"},
        contig2: {"contig_feat": "XYZ"},
    }


@pytest.mark.parametrize(
    "coordinates, is_complement, start_shift, expected",
    [
        # Case 1: No partial start or end, expect no change in non-complement strand
        # Coordinates are already correct, no need to modify anything.
        ([(11, 40)], False, 0, [(11, 40)]),
        # Case 2: Partial start, no partial end (Non-complement)
        # A shift of 1 is added to the start coordinate.
        ([(10, 40)], False, 1, [(11, 40)]),  # start_shift of 1 added to start
        # Case 2: Partial start, no partial end (Non-complement)
        # A shift of 2 is added to the start coordinate.
        ([(2, 33)], False, 2, [(4, 33)]),  # start_shift of 2 added to start
        # Case 3: No partial start, partial end (Non-complement)
        # Adjust last coordinate to make gene length a multiple of 3.
        ([(11, 41)], False, 0, [(11, 40)]),  # last end adjusted to be a multiple of 3
        # Case 3: No partial start, partial end (Non-complement)
        # Gene length is already a multiple of 3, so no changes needed.
        (
            [(11, 40)],
            False,
            0,
            [(11, 40)],
        ),  # gene length already a multiple of 3 so no change is made
        # Case 4: Partial start and end (Non-complement)
        # Both start and end need adjustment: add shift to start and adjust end to make gene length a multiple of 3.
        (
            [(10, 41)],
            False,
            1,
            [(11, 40)],
        ),  # start_shift added to start, and length adjusted
        # Case 5: Partial start and no partial end on complement strand
        # Adjust start since we are on the complement strand.
        ([(9, 40)], True, 0, [(11, 40)]),  # length adjusted
        # Case 5: No partial start but partial end on complement strand
        # Shift removed from the last end on the complement strand.
        ([(9, 40)], True, 2, [(9, 38)]),  # start_shift removed
        # Case 5: Partial start and end on complement strand
        # Adjust both start and end since we are on the complement strand, ensuring gene length is a multiple of 3.
        ([(8, 40)], True, 2, [(9, 38)]),  # start_shift removed and length adjusted
        # Case 5: Joined coordinates without partial start or end
        # Nothing to adjust as the gene is properly framed.
        ([(1, 9), (7, 12)], False, 0, [(1, 9), (7, 12)]),  # nothing to do
        # Case 5: Joined coordinates with partial start
        # Adjust the first start coordinate by the shift.
        ([(3, 9), (7, 12)], False, 1, [(4, 9), (7, 12)]),  # adjust start
        # Case 5: Joined coordinates with partial end
        # Adjust the last end to ensure the gene length is a multiple of 3.
        ([(3, 9), (7, 12)], False, 0, [(3, 9), (7, 11)]),  # adjust end
        # Case 5: Joined coordinates with partial start and partial end
        # Adjust both start and end for correct gene length and frame shift.
        ([(3, 9), (7, 12)], False, 2, [(5, 9), (7, 10)]),  # adjust start and end
        # Case 5: Joined coordinates with partial start and end on complement strand
        # Adjust both start and end on the complement strand.
        (
            [(4, 9), (7, 12)],
            True,
            2,
            [(5, 9), (7, 10)],
        ),  # adjust start and end in complement mode
        # Real tricky case from GCF_000623275.1
        ([(4681814, 4682911), (1, 1)], False, 0, [(4681814, 4682911)]),
        # ajust the end by removing one nt. In this case that remove the second part of the coordinates
        # Tricky case inspired by last case
        ([(30, 60), (1, 1)], False, 0, [(30, 59)]),
        # ajust the end by removing two nt. In this case that remove the second part of the coordinates and one nt in the first part
        # Tricky case inspired by last case
        ([(60, 60), (1, 9)], False, 1, [(1, 9)]),
        # ajust the end by removing one nt. In this case that remove the second part of the coordinates
        # Tricky case inspired by last case
        ([(60, 60), (1, 10)], False, 2, [(2, 10)]),
        # ajust the end by removing one nt. In this case that remove the second part of the coordinates
        # Very tricky case inspired by last case
        ([(59, 60), (60, 60), (1, 9)], False, 3, [(1, 9)]),
        # ajust the end by removing one nt. In this case that remove the second part of the coordinates
        # Very tricky case inspired by last case
        ([(60, 61), (1, 8)], True, 3, [(61, 61), (1, 5)]),  #
    ],
)
def test_fix_partial_gene_coordinates(
    coordinates, is_complement, start_shift, expected
):
    result = fix_partial_gene_coordinates(coordinates, is_complement, start_shift)
    assert result == expected


def test_fix_partial_gene_coordinates_with_wrong_coordinates():
    with pytest.raises(ValueError):
        fix_partial_gene_coordinates(
            coordinates=[(1, 1)], is_complement=False, start_shift=0
        )  # gene is too small, the length adjustement at the end lead to no gene

    with pytest.raises(ValueError):
        fix_partial_gene_coordinates(
            [(1, 1)], False, 1
        )  # gene is too small, the length adjustement at the start lead to no gene

    with pytest.raises(ValueError):
        fix_partial_gene_coordinates(
            [(60, 60), (1, 1)], False, 1
        )  # gene is too small, the length adjustement at the start and at the end lead to no gene

    with pytest.raises(ValueError):
        fix_partial_gene_coordinates([], False, 1)  # chevron in inner position


@pytest.mark.parametrize(
    "coordinates, shift, expected",
    [
        ([(11, 40)], 0, [(11, 40)]),
        ([(1, 2)], 1, [(2, 2)]),
        ([(1, 1), (1, 4)], 1, [(1, 4)]),
        ([(1, 1), (1, 1), (1, 4)], 2, [(1, 4)]),
        ([(1, 1), (1, 2), (1, 4)], 2, [(2, 2), (1, 4)]),
    ],
)
def test_shift_start_coordinates(coordinates, shift, expected):
    result = shift_start_coordinates(coordinates, shift)
    assert result == expected


@pytest.mark.parametrize(
    "coordinates, shift, expected",
    [
        ([(11, 40)], 0, [(11, 40)]),
        ([(1, 2)], 1, [(1, 1)]),
        ([(1, 1), (1, 4)], 1, [(1, 1), (1, 3)]),
        ([(1, 4), (4, 4), (4, 4)], 2, [(1, 4)]),
        ([(18, 18), (1, 4)], 4, [(18, 18)]),
    ],
)
def test_shift_end_coordinates(coordinates, shift, expected):
    result = shift_end_coordinates(coordinates, shift)
    assert result == expected


def test_check_sequence_tuple_valid():
    name, sequence = check_sequence_tuple("seq1", "ATGC")
    assert name == "seq1"
    assert sequence == "ATGC"


def test_check_sequence_tuple_empty_name():
    with pytest.raises(ValueError):
        check_sequence_tuple("", "ATGC")


def test_check_sequence_tuple_empty_sequence():
    with pytest.raises(ValueError):
        check_sequence_tuple("seq1", "")


def test_parse_fasta_valid():
    fasta_data = ">seq1\nATGC\n>seq2\nGCTA"

    result = list(parse_fasta(fasta_data.split("\n")))

    assert result == [("seq1", "ATGC"), ("seq2", "GCTA")]


def test_parse_fasta_empty_sequence():
    fasta_data = ">seq1\n>seq2\nGCTA"
    with pytest.raises(ValueError):
        list(parse_fasta(fasta_data.split("\n")))


def test_parse_fasta_no_header():
    fasta_data = "seq1\nATGC\nseq2\nGCTA".split("\n")
    with pytest.raises(ValueError):
        list(parse_fasta(fasta_data))
