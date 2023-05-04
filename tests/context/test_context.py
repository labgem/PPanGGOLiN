import pytest
from ppanggolin.context.searchGeneContext import extract_contig_window


def test_extract_contig_window():
    assert extract_contig_window(contig_length=15, positions_of_interest={8}, window_size=1) == [(7,9)]

    # check that extracted window is inside contig limit
    assert extract_contig_window(contig_length=16, positions_of_interest={15}, window_size=4) == [(11,15)]

    assert extract_contig_window(contig_length=10, positions_of_interest={2, 8}, window_size=2) == [(0,4), (6,9)]
    
    # 12 window is (9,15)
    # 19  window is (16,22)
    # so when 12 and 19 are of interest window merge (9,22)
    assert extract_contig_window(contig_length=200, positions_of_interest={12}, window_size=3) == [(9,15)]
    assert extract_contig_window(contig_length=200, positions_of_interest={19}, window_size=3) == [(16,22)]
    assert extract_contig_window(contig_length=200, positions_of_interest={12, 19}, window_size=3) == [(9,22)]

    assert extract_contig_window(contig_length=10, positions_of_interest={2, 5, 8}, window_size=2) == [(0,9)]

def test_extract_contig_window_with_circular_contig():
    # # check that circularity is properly taken into account
    assert extract_contig_window(contig_length=12, positions_of_interest={1}, window_size=2, is_circular=True) == [(0,3), (11,11)]
    assert extract_contig_window(contig_length=12, positions_of_interest={1}, window_size=3, is_circular=True) == [(0,4), (10,11)]
    assert extract_contig_window(contig_length=12, positions_of_interest={10}, window_size=3, is_circular=True) == [(0,1), (7,11)]

    assert list(extract_contig_window(contig_length=12, positions_of_interest={6}, window_size=6, is_circular=True)) == [(0,11)]
    assert list(extract_contig_window(contig_length=12, positions_of_interest={1}, window_size=6, is_circular=True)) == [(0,11)]
    assert list(extract_contig_window(contig_length=12, positions_of_interest={1}, window_size=6, is_circular=False)) == [(0,7)]

    assert list(extract_contig_window(contig_length=12, positions_of_interest={0, 9}, window_size=2, is_circular=False)) == [(0,2), (7,11)]

def test_extract_contig_window_out_of_range():
    with pytest.raises(IndexError):
        extract_contig_window(contig_length=15, positions_of_interest={15}, window_size=1)

    with pytest.raises(IndexError):
        extract_contig_window(contig_length=15, positions_of_interest={-1}, window_size=1)