#! /usr/bin/env python3

import pytest

from genome import RNA

"""
"""
@pytest.fixture()
def o_rna():
    return RNA(4)

def test_cstr():
    ID  = 4
    o_rna = RNA(ID)
    assert isinstance(o_rna, RNA)
    for attr in "ID", "is_fragment", "type":
        assert hasattr(o_rna, attr)
    assert o_rna.ID == ID
    assert o_rna.is_fragment == False
    assert o_rna.type == ""

def test_fill_annotations(o_rna):
    start,stop = 1,9
    strand = "plus"
    o_rna.fill_annotations(start, stop, strand)
    for attr in 'start', 'stop', 'strand', \
                'type', 'product', 'name':
        assert hasattr(o_rna, attr)
    assert o_rna.start == start
    assert o_rna.stop  == stop
    assert o_rna.strand== strand
    assert o_rna.type  == ''
    assert o_rna.name  == ''
    assert o_rna.product == ''

    gene_type = "inconnu"
    name      = "Eugène"
    product   = "va savoir"
    o_rna.fill_annotations(start, stop, strand, gene_type, name, product)
    assert o_rna.type == gene_type
    assert o_rna.name == name
    assert o_rna.product == product

    # what if start or stop < 0 ?
    #   stop < start
    # start/stop cannot int() ?
    # position not int

def test_fill_parents(o_rna):
    org = "toto"
    ctg = 99
    o_rna.fill_parents(org, ctg)
    for attr in 'organism', 'contig':
        assert hasattr(o_rna, attr)
    assert o_rna.organism == org
    assert o_rna.contig == ctg

def test_add_dna(o_rna):
    dna = "test adn"
    o_rna.add_dna(dna)
    assert hasattr(o_rna, 'dna')
    o_rna.dna = dna

    dna = 123
    with pytest.raises(TypeError):
        o_rna.add_dna(dna)
