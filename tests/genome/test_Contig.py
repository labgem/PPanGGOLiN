#! /usr/bin/env python3
import random

import pytest

from ppanggolin.genome import Contig, Gene, RNA


@pytest.fixture()
def o_ctg():
    return Contig("toto")


def test_cstr():
    name = 4
    o_ctg = Contig(name)
    assert isinstance(o_ctg, Contig)
    for attr in "name", "is_circular", "RNAs":
        assert hasattr(o_ctg, attr)
    assert o_ctg.name == name
    assert o_ctg.is_circular is False
    assert o_ctg.RNAs == set()

    o_ctg = Contig(name, True)
    assert o_ctg.is_circular is True


def test_str():
    name = "ppoiu"
    o_ctg = Contig(name)
    assert str(o_ctg) == name


def test_add_rna(o_ctg):
    with pytest.raises(TypeError):
        o_ctg.add_rna(33)

    l_rnas = []
    for i in "abdc":
        o_rna = RNA(i)
        o_ctg.add_rna(o_rna)
        l_rnas.append(o_rna)
    assert o_ctg.RNAs == set(l_rnas)


@pytest.fixture()
def l_genes():
    l_genes = []
    for i in range(6, -1, -1):  # Create 7 Gene
        o_gene = Gene(i)
        o_gene.fill_annotations(start=i*10, stop=i*10 - 1, strand='+', position=i)
        l_genes.append(o_gene)

    return l_genes


def test_add_gene(o_ctg, l_genes):
    with pytest.raises(TypeError):
        o_ctg.add_gene(33)

    # gene must have a position before beeing added.
    with pytest.raises(TypeError):
        o_ctg.add_gene(Gene(33))

    for o_gene in l_genes:
        o_ctg.add_gene(o_gene)

    assert o_ctg.genes == sorted(l_genes, key=lambda x: x.position)


def test_iterator_behavior(o_ctg, l_genes):
    # FIXME: is there a better way to check this ?
    assert iter(o_ctg)

    for o_gene in l_genes:
        o_ctg.add_gene(o_gene)

    l_ = [o_gene for o_gene in o_ctg]
    assert l_ == sorted(l_genes, key=lambda x: x.start)
