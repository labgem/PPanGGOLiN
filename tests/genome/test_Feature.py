#! /usr/bin/env python3

import pytest

from ppanggolin.genome import Feature


def test_cstr():
    identifier = 4
    o_feature = Feature(identifier)
    assert isinstance(o_feature, Feature)
    for attr in "ID", "is_fragment", "type":
        assert hasattr(o_feature, attr)
    assert o_feature.ID == identifier
    assert o_feature.is_fragment is False
    assert o_feature.type == ""


@pytest.fixture()
def o_feature():
    return Feature(4)


def test_fill_annotations(o_feature):
    start, stop = 1, 9
    strand = "plus"
    o_feature.fill_annotations(start, stop, strand)
    for attr in 'start', 'stop', 'strand', \
                'type', 'product', 'name':
        assert hasattr(o_feature, attr)
    assert o_feature.start == start
    assert o_feature.stop == stop
    assert o_feature.strand == strand
    assert o_feature.type == ''
    assert o_feature.name == ''
    assert o_feature.product == ''

    gene_type = "inconnu"
    name = "Eugène"
    product = "va savoir"
    o_feature.fill_annotations(start, stop, strand, gene_type, name, product)
    assert o_feature.type == gene_type
    assert o_feature.name == name
    assert o_feature.product == product

    # what if start or stop < 0 ?
    #   stop < start
    # start/stop cannot int() ?
    # position not int


def test_fill_parents(o_feature):
    org = "toto"
    ctg = 99
    o_feature.fill_parents(org, ctg)
    for attr in 'organism', 'contig':
        assert hasattr(o_feature, attr)
    assert o_feature.organism == org
    assert o_feature.contig == ctg


def test_add_dna(o_feature):
    dna = "test adn"
    o_feature.add_dna(dna)
    assert hasattr(o_feature, 'dna')
    o_feature.dna = dna

    dna = 123
    with pytest.raises(TypeError):
        o_feature.add_dna(dna)
