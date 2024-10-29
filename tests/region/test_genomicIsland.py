from ppanggolin.utils import (
    find_consecutive_sequences,
    find_region_border_position,
    get_consecutive_region_positions,
)
import pytest


def test_find_consecutive_sequences_single_sequence():
    sequence = [1, 2, 3, 4, 5]
    assert find_consecutive_sequences(sequence) == [[1, 2, 3, 4, 5]]


def test_find_consecutive_sequences_multiple_sequences():
    sequence = [1, 2, 3, 7, 8, 9, 11, 12, 13, 0]
    assert find_consecutive_sequences(sequence) == [
        [0, 1, 2, 3],
        [7, 8, 9],
        [11, 12, 13],
    ]


def test_find_region_border_position_single_sequence():
    region_positions = [1, 2, 3, 4, 5]
    contig_length = 10
    assert find_region_border_position(region_positions, contig_length) == (1, 5)


def test_find_region_border_position_edge_overlap():
    region_positions = [0, 1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 14]
    contig_length = 15
    assert find_region_border_position(region_positions, contig_length) == (7, 3)


def test_find_region_border_position_empty_sequence():
    region_positions = []
    contig_length = 10
    with pytest.raises(ValueError):
        find_region_border_position(region_positions, contig_length)


def test_find_region_border_position_multiple_fragments():
    region_positions = [1, 3, 4, 5, 8, 9]
    contig_length = 10
    with pytest.raises(ValueError):
        find_region_border_position(region_positions, contig_length)


def test_find_region_border_position_fragmented_but_no_zero():
    # region is in two pieces but it miss position 0 to be correct overlap
    region_positions = [8, 9, 1, 2, 3, 4]
    contig_length = 10
    with pytest.raises(ValueError):
        find_region_border_position(region_positions, contig_length)


def test_find_region_border_position_fragmented_but_wrong_max_po():
    # region does not reach the end of the contig. It misses position 9
    region_positions = [7, 8, 0, 1, 2, 3, 4]
    contig_length = 10
    with pytest.raises(ValueError):
        find_region_border_position(region_positions, contig_length)


def test_get_consecutive_region_positions_regular():
    region_positions = [2, 3, 4, 5, 6]
    contig_length = 15
    assert get_consecutive_region_positions(region_positions, contig_length) == [
        [2, 3, 4, 5, 6]
    ]


def test_get_consecutive_region_positions_overlap():
    region_positions = [0, 1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 14]
    contig_length = 15
    assert get_consecutive_region_positions(region_positions, contig_length) == [
        [7, 8, 9, 10, 11, 12, 13, 14],
        [0, 1, 2, 3],
    ]


def test_get_consecutive_region_positions_all_genes():
    region_positions = [4, 5, 6, 7, 0, 1, 2, 3]
    contig_length = 8
    assert get_consecutive_region_positions(region_positions, contig_length) == [
        [0, 1, 2, 3, 4, 5, 6, 7]
    ]
