#! /usr/bin/env python3

import pytest
from collections import defaultdict

from ppanggolin.genome import Gene
from ppanggolin.pangenome import Edge, GeneFamily

"""
"""
def test_cstr_error():
    o_src = Gene('source')
    o_tgt = Gene('target')
    # genes should have a family
    with pytest.raises(Exception):
        o_edge = Edge(o_src, o_tgt)

    o_family = GeneFamily(None, None)
    o_family.addGene(o_src)
    # both genes sould have a family
    with pytest.raises(Exception):
        o_edge = Edge(o_src, o_tgt)

    # gene should belong to the same organism
    o_family.addGene(o_tgt)
    o_src.fill_parents("",None)
    o_tgt.fill_parents(None,None)
    with pytest.raises(Exception):
        o_edge = Edge(o_src, o_tgt)


def test_cstr():
    o_src = Gene('source')
    o_tgt = Gene('target')

    # set organism and contig to None.
    o_src.fill_parents(None,None)
    o_tgt.fill_parents(None,None)

    # define the None GeneFamily, and add the 2 genes to it.
    o_family = GeneFamily(None, None)
    o_family.addGene(o_src)
    o_family.addGene(o_tgt)

    o_edge = Edge(o_src, o_tgt)
    assert isinstance(o_edge, Edge)

    assert o_edge.source == o_src.family
    assert o_edge.target == o_tgt.family
    assert dict(o_edge.organisms) == { None: [(o_src, o_tgt)] }


@pytest.fixture()
def make_gene_pair():
    def _make_gene_pair(org, gene_id1, gene_id2):
        """create 2 genes from org.
            each gene belong to its own family."""
        lo_genes = []
        for k in gene_id1, gene_id2:
            o_gene = Gene(k)
            o_gene.fill_parents(org,None)

            lo_genes.append(o_gene)

            o_family = GeneFamily(k,k)
            o_family.addGene(o_gene)

        return tuple(lo_genes)

    return _make_gene_pair


@pytest.fixture()
def o_edge(make_gene_pair):
    p = make_gene_pair("org", "src", "tgt")
    return Edge(*p)


def test_addGenes(make_gene_pair):
    p1 = make_gene_pair("org1", "s1", "t1")
    p2 = make_gene_pair("org1", "s2", "t1")
    p3 = make_gene_pair("org2", "s1", "t2")
    p4 = make_gene_pair("org2", "s1", "s2")
    # org1: s1,s2 -- t1
    # org2: s1 -- t2,s2

    o_edge = Edge(*p1)
    o_edge.addGenes(*p2)
    o_edge.addGenes(*p3)
    o_edge.addGenes(*p4)
    assert set(o_edge.organisms.keys()) == set(["org1", "org2"])
    assert o_edge.organisms["org1"] == [p1,p2]
    assert o_edge.organisms["org2"] == [p3,p4]


@pytest.fixture()
def filled_edge(make_gene_pair):
    # Note that the same edge here links 4 families.
    p1 = make_gene_pair("org1", "s1", "t1")
    p2 = make_gene_pair("org1", "s2", "t1")
    p3 = make_gene_pair("org2", "s1", "t2")
    p4 = make_gene_pair("org2", "s1", "s2")
    # org1: s1,s2 -- t1
    # org2: s1 -- t2,s2

    o_edge = Edge(*p1)
    o_edge.addGenes(*p2)
    o_edge.addGenes(*p3)
    o_edge.addGenes(*p4)

    return o_edge


def test_getOrgDict(o_edge, filled_edge):
    assert o_edge.getOrgDict() == o_edge.organisms
    assert filled_edge.getOrgDict() == filled_edge.organisms


def test_genePairs(make_gene_pair):
    # cannot use filled_edge because i need access to pairs.
    p1 = make_gene_pair("org1", "s1", "t1")
    p2 = make_gene_pair("org1", "s2", "t1")
    p3 = make_gene_pair("org2", "s1", "t2")
    p4 = make_gene_pair("org2", "s1", "s2")
    # org1: s1,s2 -- t1
    # org2: s1 -- t2,s2

    o_edge = Edge(*p1)
    o_edge.addGenes(*p2)
    o_edge.addGenes(*p3)
    o_edge.addGenes(*p4)

    # 'set' because the order is not guaranted due to '.values()'.
    l_pairs = o_edge.genePairs
    assert set(l_pairs) == set([p1,p2,p3,p4])