from ppanggolin.formats.writeFlatGenomes import convert_overlapping_coordinates_for_gff


def test_convert_overlapping_coordinates_for_gff():
    # test case where coordinates are have no frameshift and no edge overlap
    coordinates = [(7, 10)]
    contig_length = 15
    assert convert_overlapping_coordinates_for_gff(coordinates, contig_length) == [
        (7, 10)
    ]

    # test case where coordinates are simply a frameshift with not edge overlap
    coordinates = [(1, 5), (7, 10)]
    contig_length = 15
    assert convert_overlapping_coordinates_for_gff(coordinates, contig_length) == [
        (1, 5),
        (7, 10),
    ]

    # Test case where the gene overlaps. Coordinates are adjusted for GFF format
    coordinates = [(4, 8), (1, 2)]
    contig_length = 8
    assert convert_overlapping_coordinates_for_gff(coordinates, contig_length) == [
        (4, 10)
    ]

    # Test case where coordinates overlap and has frameshift
    coordinates = [(4, 8), (7, 13), (1, 2)]
    contig_length = 13
    assert convert_overlapping_coordinates_for_gff(coordinates, contig_length) == [
        (4, 8),
        (7, 15),
    ]

    # Test case where coordinates overlap and has frameshift at the edge of the contig
    coordinates = [(12, 18), (1, 4)]
    contig_length = 20
    assert convert_overlapping_coordinates_for_gff(coordinates, contig_length) == [
        (12, 18),
        (21, 24),
    ]

    # Test case where coordinates overlap and has frameshift at the edge of the contig
    coordinates = [(12, 20), (2, 4)]
    contig_length = 20
    assert convert_overlapping_coordinates_for_gff(coordinates, contig_length) == [
        (12, 20),
        (22, 24),
    ]

    # Test case where coordinates have just the last nt of the contig and the rest is at th ebegining of the contig
    coordinates = [(20, 20), (1, 10)]
    contig_length = 20
    assert convert_overlapping_coordinates_for_gff(coordinates, contig_length) == [
        (20, 30)
    ]
