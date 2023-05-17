#! /usr/bin/env python3

import pytest

from ppanggolin.RGP import rgp_cluster
from ppanggolin.region import Region

from test_Region import l_genes, o_org, o_contig


@pytest.fixture
def identical_rgps(l_genes):
    rgp1 = Region("one")
    rgp2 = Region('two')
    rgp3 = Region('three')
    
    # the three rgp have the gene content. 
    # in terms of family they are identical
    for g in l_genes:
        rgp1.append(g)
        rgp2.append(g)
        rgp3.append(g)
    
    return [rgp1, rgp2, rgp3]

@pytest.fixture
def RGP_a(l_genes):
    rgp = Region("A")
    
    for g in l_genes[:8]:
        rgp.append(g)

    return rgp


def o_region():
    return Region(4)

def test_compute_grr():
    set1 = {1, 2, 3, 4, 5}
    set2 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

    assert rgp_cluster.compute_grr(set1, set2, min) == 1.0
    assert rgp_cluster.compute_grr(set1, set2, max) == 0.5

    
def test_dereplicate_rgp(identical_rgps):

    rgp1 = identical_rgps[0]
    assert rgp_cluster.dereplicate_rgp([rgp1]) == {rgp1:{rgp1}}

    assert rgp_cluster.dereplicate_rgp(identical_rgps) == {identical_rgps[0]:set(identical_rgps)}


def test_compute_rgp_metric():
    
    rgp_pair = (1, 2)
    rgp_to_families = {1: {"A", 'B'},
                       2: {"A", "B", 'C'} }
    rgp_to_contigborder = {1: True,
                           2: False }
    # min_grr
    min_result = rgp_cluster.compute_rgp_metric(rgp_pair, rgp_to_families, rgp_to_contigborder, 0.8, "min_grr")
    expected_min_grr = (1,2, {'incomplete_aware_grr':2/2, 
                            "min_grr":2/2, 
                            'max_grr':2/3,
                            'shared_family':2})
    assert min_result == expected_min_grr

    # incomplete_aware_grr: same as min grr as rgp1 is incomplete 
    incomplete_aware_result = rgp_cluster.compute_rgp_metric(rgp_pair, rgp_to_families, rgp_to_contigborder, 0.8, "incomplete_aware_grr")

    assert incomplete_aware_result == expected_min_grr


    # max grr is below cutoff so None is returned
    assert rgp_cluster.compute_rgp_metric(rgp_pair, rgp_to_families, rgp_to_contigborder, 0.8, "max_grr") == None


def test_get_rgp_info_dict(RGP_a):
    rgp_id = RGP_a.ID
    region_to_spot={RGP_a:6}
    region_info = {"contig": 1,
                    'organism': "toto",
                    "name": "A",
                    "genes_count": 8,
                    "is_contig_border": True,
                    "is_whole_contig": True,
                    "spot_id": '6',
                    "families_count":8
                    }
    
    assert rgp_cluster.get_rgp_info_dict([RGP_a], region_to_spot) == {rgp_id:region_info}


