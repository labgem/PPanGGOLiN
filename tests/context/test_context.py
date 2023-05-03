import pytest
from ppanggolin.context.searchGeneContext import extract_contig_window


def test_extract_contig_window():
    assert list(extract_contig_window(contig_length=15, positions_of_interest={8}, window_size=1)) == [(7,9)]

    # check that extracted window is inside contig limit
    assert list(extract_contig_window(contig_length=16, positions_of_interest={15}, window_size=4)) == [(11,15)]

    assert list(extract_contig_window(contig_length=10, positions_of_interest={2, 8}, window_size=2)) == [(0,4), (6,9)]

    assert list(extract_contig_window(contig_length=10, positions_of_interest={2, 5, 8}, window_size=2)) == [(0,9)]

    # # check that circularity is properly taken into account
    # assert list(extract_contig_window(contig_length=11, positions_of_interest={0}, window_size=2, is_circular=True)) == [(9,2)]

    # assert list(extract_contig_window(contig_length=11, positions_of_interest={0, 9}, window_size=2, is_circular=True)) == [(7,2)]

    # assert list(extract_contig_window(contig_length=11, positions_of_interest={0, 9}, window_size=2, is_circular=False)) == [(0,2), (7,11)]


    with pytest.raises(IndexError):
        list(extract_contig_window(contig_length=15, positions_of_interest={15}, window_size=1))

    with pytest.raises(IndexError):
        list(extract_contig_window(contig_length=15, positions_of_interest={-1}, window_size=1))