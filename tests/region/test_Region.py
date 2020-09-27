#! /usr/bin/env python3

import pytest
from random import choices, randint, sample

from ppanggolin.region import Region
from ppanggolin.genome import Gene, Contig, Organism, RNA

def test_cstr():
    ID = 4
    o_region = Region(ID)
    assert isinstance(o_region, Region)
    for attr in "genes", "name", "score":
        assert hasattr(o_region, attr)
    
    assert o_region.score == 0
    assert o_region.name == ID
    assert o_region.genes == []

@pytest.fixture
def o_region():
    return Region(4)

@pytest.fixture
def o_org():
    return Organism("toto")

@pytest.fixture
def o_contig():
    return Contig(1)

@pytest.fixture
def o_rna(o_contig):
    o_rna = RNA("Ah")
    o_rna.fill_annotations(35,45, "-")
    o_contig.addRNA(o_rna)
    return o_rna

@pytest.fixture
def l_genes(o_org, o_contig):
    """ creates a small testing context, with 4 CDS, 1 RNA that are all on the same contig in the same organism"""
    l_genes = []
    c=10
    for i, gene_id in enumerate(["toto","tata","titi","tutu"]):
        gene = Gene(gene_id)
        gene.fill_annotations(c,c+30, "+",position=i)
        gene.fill_parents(o_org, o_contig)
        o_contig.addGene(gene)
        gene.family = gene_id
        l_genes.append(gene)
        c+=35
    return l_genes

def test_append(l_genes, o_region):
    for gene in l_genes:
        o_region.append(gene)
    
    assert set(o_region.genes) == set(l_genes)

def test_append__error(o_region):
    """append should raise a TypeError is used with non Gene param."""
    with pytest.raises(TypeError):
        o_region.append(42)

def test_properties(l_genes, o_region, o_org, o_contig):
    s_families = set()
    for gene in l_genes:
        o_region.append(gene)
        s_families.add(gene.family)

    #checking properties sanity
    assert o_region.start == o_region.startGene.start
    assert o_region.stop == o_region.stopGene.stop
    assert o_region.isWholeContig is True
    assert o_region.isContigBorder is True
    assert o_region.contig == o_contig
    assert o_region.organism == o_org
    assert o_region.families == s_families

def test_getRNAs(o_rna, o_region, l_genes):
    for gene in l_genes:
        o_region.append(gene)
    assert set(o_region.getRNAs()) == set([o_rna])

def test_hash(o_region):
    """ a hash function returns an integer"""
    # the same int if called twice on the same object
    h = hash(o_region)
    assert isinstance(h, int)
    assert h == hash(o_region)

    # different ints if called on objects representing the same entity
    name = "charming"
    assert hash(Region(name)) != hash(Region(name))

def test_equality(o_region, l_genes):
    """2 regions are equals if they contains the same list of genes."""
    for gene in l_genes:
        o_region.append(gene)

    # not the same list => False
    o_other = Region("other")
    assert o_region != o_other

    # the exact same list => True
    o_other = Region("other")
    for gene in l_genes:
        o_other.append(gene)
    assert o_region == o_other

    # the same list in reverse order => True
    o_other = Region("other")
    for gene in reversed(l_genes):
        o_other.append(gene)
    assert o_region == o_other

def test_equality__error(o_region):
    """equality raises error if not compared to another Region"""
    with pytest.raises( TypeError ):
        o_region == 42
