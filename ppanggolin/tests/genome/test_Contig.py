#! /usr/bin/env python3

import pytest

from genome import Contig, Gene, RNA

"""
"""
@pytest.fixture()
def o_ctg():
    return Contig("toto")

def test_cstr():
    name  = 4
    o_ctg = Contig(name)
    assert isinstance(o_ctg, Contig)
    for attr in "name", "is_circular", "RNAs":
        assert hasattr(o_ctg, attr)
    assert o_ctg.name == name
    assert o_ctg.is_circular == False
    assert o_ctg.RNAs == set()

    o_ctg = Contig(name, True)
    assert o_ctg.is_circular == True

def test_str():
    name = "ppoiu"
    o_ctg = Contig(name)
    assert str(o_ctg) == name

def test_addSequence(o_ctg):
    seq = "sequence"
    o_ctg.addSequence(seq)
    assert hasattr(o_ctg, seq)
    assert o_ctg.sequence == seq

    with pytest.raises(TypeError):
        o_ctg.addSequence(33)

def test_addRNA(o_ctg):
    with pytest.raises(TypeError):
        o_ctg.addRNA(33)

    l_rnas = []
    for i in "abdc":
        o_rna = RNA(i)
        o_ctg.addRNA(o_rna)
        l_rnas.append(o_rna)
    assert o_ctg.RNAs == set(l_rnas)

@pytest.fixture()
def l_genes():
    l_genes = []
    for i in range(6,0,-1):
        o_gene = Gene(i)
        o_gene.fill_annotations(i,i,i,position=i-1)
        l_genes.append(o_gene)

    return l_genes

def test_addGene(o_ctg, l_genes):
    with pytest.raises(TypeError):
        o_ctg.addGene(33)

    # gene must have a position before beeing added.
    with pytest.raises(TypeError):
        o_ctg.addGene(Gene(33))

    for o_gene in l_genes:
        o_ctg.addGene(o_gene)

    assert o_ctg.genes == sorted(l_genes, key=lambda x:x.position)

def test_iterator_behavior(o_ctg, l_genes):
    # FIXME: is there a better way to check this ?
    assert iter(o_ctg)

    for o_gene in l_genes:
        o_ctg.addGene(o_gene)

    l_ = [ o_gene for o_gene in o_ctg ]
    assert l_ == sorted(l_genes, key=lambda x:x.start)
