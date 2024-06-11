import pytest
from pathlib import Path
from ppanggolin.annotate.annotate import extract_positions, read_anno_file, parse_contig_header_lines, parse_gbff_by_contig, parse_feature_lines, parse_dna_seq_lines, read_org_gbff, combine_contigs_metadata
    


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
    script_path = Path(__file__).resolve()
    ppanggolin_main_dir = script_path.parent.parent.parent

    genome_path = ppanggolin_main_dir / "testingDataset/GBFF/GCF_000026905.1_ASM2690v1_genomic.gbff.gz"
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

    genome_path = ppanggolin_main_dir / "testingDataset/GBFF/GCF_002776845.1_ASM277684v1_genomic.gbff.gz"
    circular_contigs = []
    genome_name = "GCF_002776845" 
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


def test_with_joined_genes(genome_data_with_joined_genes):
    genome_name, genome_path, circular_contigs = genome_data_with_joined_genes
    use_pseudogene = True
    genome, _ = read_anno_file(genome_name, genome_path, circular_contigs, use_pseudogene)

    # this genome has 2 genes that are joined.
    assert genome.number_of_genes() == 917 
    
def test_read_org_gbff(genome_data_with_joined_genes):
    genome_name, genome_path, circular_contigs = genome_data_with_joined_genes
    genome, _ = read_org_gbff(genome_name, genome_path, circular_contigs, pseudo=True)

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

    assert parsed_header ==  {
            "LOCUS"     :  "NC_022109            1041595 bp    DNA     circular CON 24-MAR-2017",
            "DEFINITION" : "Chlamydia trachomatis strain D/14-96 genome.",
            "VERSION"   :  "NC_022109.1",
            "SOURCE"     : "Chlamydia trachomatis",
            "ORGANISM" : "Chlamydia trachomatis\nBacteria; Chlamydiae; Chlamydiales; Chlamydiaceae;\nChlamydia/Chlamydophila group; Chlamydia."
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
    assert list(feature_2) == [{"type":"source", "location":"1..2041595", "organism":'Chlamydia trachomatis', "mol_type":"genomic DNA"}]

    assert sequence_2 == "AAACCGGGTTCCAAATTTGGGGCCCCTTTT"



# Define test cases
@pytest.mark.parametrize("input_lines, expected_output", [
    (["     gene            123..456",
      "                     /locus_tag=ABC123",
      "                     /note=Some note",
      "     CDS             789..1011",
      "                     /protein_id=DEF456",
      "                     /translation=ATGCTAGCATCG"], 
     [{"type": "gene",
       "location": "123..456",
       "locus_tag": "ABC123",
       "note": "Some note"},
      {"type": "CDS","location": "789..1011", "protein_id": "DEF456", "translation": "ATGCTAGCATCG"}]),
    (["     gene            123..456", 
      "                     /locus_tag=ABC123", 
      "                     /note=Some note", 
      "     CDS             789..1011",
      '                     /protein_id="DEF456"', 
      '                     /translation="ATGCTAGCATCG"', 
      "     gene            789..1011",
      "                     /locus_tag=DEF789",
      "                     /note=Another note"], 
     [{"type": "gene", "location": "123..456", "locus_tag": "ABC123", "note": "Some note"},
      {"type": "CDS", "location": "789..1011", "protein_id": "DEF456", "translation": "ATGCTAGCATCG"},
      {"type": "gene", "location": "789..1011", "locus_tag": "DEF789", "note": "Another note"}]),
    # Add more test cases as needed
])
def test_parse_feature_lines(input_lines, expected_output):
    print(expected_output)
    assert list(parse_feature_lines(input_lines)) == expected_output


def test_parse_dna_seq_lines():
    lines =["        1 aaacc gggtt",
            "       11 ccaaa tttgg",
            "       21 ggccc ctttt"]
    
    assert parse_dna_seq_lines(lines) == "AAACCGGGTTCCAAATTTGGGGCCCCTTTT"


def test_combine_contigs_metadata():

    contig_to_metadata = {
        "contig1":{"sp":"spA", "strain":"123", "contig_feat":"ABC"},
        "contig2":{"sp":"spA", "strain":"123", "contig_feat":"XYZ"},
        "contig3":{"sp":"spA", "strain":"123"}
    }

    genome_metadata, contig_metadata = combine_contigs_metadata(contig_to_metadata)

    assert genome_metadata == {"sp":"spA", "strain":"123"}
    assert contig_metadata == {"contig1":{"contig_feat":"ABC"}, "contig2":{"contig_feat":"XYZ"},}
