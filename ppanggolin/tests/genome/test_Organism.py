#! /usr/bin/env python3

import pytest
from random import randint

from ppanggolin.genome import Contig, Gene, Organism

"""
"""
def test_cstr():
    name  = 4
    o_org = Organism(name)
    assert isinstance(o_org, Organism)
    assert hasattr(o_org, "name")
    assert o_org.name == name

def test_str():
    name = "ppoiu"
    o_org = Organism(name)
    assert str(o_org) == name

@pytest.fixture()
def o_org():
    return Organism("toto")

def test_addContig(o_org):
    #FIXME: shouldn't the method be called getContig ?
    o_ctg = o_org.addContig('i')
    assert isinstance(o_ctg, Contig)

@pytest.fixture()
def t_filled_org(o_org):
    from random import randint
    n = 0
    for k in "azerty'":
        o_ctg = o_org.addContig(k)
        for i in range(randint(0,5)):
            o_gene = Gene(k+"-"+str(i))
            o_gene.fill_annotations(6,1,k,position=i)
            o_ctg.addGene(o_gene)
            n += 1

    return o_org,n

def test_families(t_filled_org):
    o_filled_org, n = t_filled_org

    # families is never set
    assert o_filled_org.families == {None}

def test_number_of_genes(t_filled_org):
    o_filled_org, n = t_filled_org

    assert o_filled_org.number_of_genes() == n

def get_genes():
    for i in range(randint(0,5)):
        o_gene = Gene(str(i))
        start  = randint(0,100)
        stop   = randint(0,100)
        o_gene.fill_annotations(start, stop,'x',position=i)
        yield o_gene

def test_contigs(o_org):
    l_contigs= []
    for k in "azer'":
        o_ctg  = o_org.addContig(k)
        for o_gene in get_genes():
            o_ctg.addGene(o_gene)
        l_contigs.append(o_ctg)

    assert list(o_org.contigs) == l_contigs

def test_genes(o_org):
    o_ctg  = o_org.addContig("scrap")
    for o_gene in get_genes():
        o_ctg.addGene(o_gene)

    assert list(o_org.genes) == o_ctg.genes

    # FIXME: find a way to test when several contigs.
    #   => order of contig is not predictable.
