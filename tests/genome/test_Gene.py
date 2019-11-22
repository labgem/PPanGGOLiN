#! /usr/bin/env python3

import pytest

from ppanggolin.genome import Feature, Gene

def test_cstr():
    """ By checking o_gene is a Feature, I rely on Feature tests."""
    ID  = 4
    o_gene = Gene(ID)
    assert isinstance(o_gene, Feature)
    assert isinstance(o_gene, Gene)

    for attr in "position", "family":
        assert hasattr(o_gene, attr)
    assert o_gene.position is None
    assert o_gene.family is None

def test_str():
    ID  = "un truc"
    o_gene = Gene(ID)
    assert str(o_gene) == ID

@pytest.fixture()
def o_gene():
    return Gene(4)

def test_fill_annotations_defaults(o_gene):
    start,stop = 1,9
    strand = "plus"
    o_gene.fill_annotations(start, stop, strand)
    for attr in "position", "genetic_code":
        assert hasattr(o_gene, attr)

    assert o_gene.position is None
    assert o_gene.genetic_code == 11

def test_fill_annotations(o_gene):
    start,stop = 1,9
    strand = "plus"
    position = "44"
    genetic_code = "le code"
    o_gene.fill_annotations(start, stop, strand, \
                                position=position,
                                genetic_code=genetic_code)
    assert o_gene.position == position
    assert o_gene.genetic_code == genetic_code

def test_add_protein_error(o_gene):
    with pytest.raises(TypeError):
        o_gene.add_protein(42)

def test_add_protein(o_gene):
    prot = "une jolie protéïne, même avec des caractères bizarres ;)"
    o_gene.add_protein(prot)
    assert o_gene.protein == prot
