#! /usr/bin/env python3

import pytest
from ppanggolin.context.searchGeneContext import (
    extract_contig_window,
    get_n_next_genes_index,
    add_edges_to_context_graph,
    compute_gene_context_graph,
)

from ppanggolin.geneFamily import GeneFamily
from ppanggolin.genome import Gene, Contig, Organism

import networkx as nx


def test_extract_contig_window():
    # TODO try to use @pytest.mark.parametrize to test different combinations
    assert extract_contig_window(
        contig_size=15, positions_of_interest={8}, window_size=1
    ) == [(7, 9)]

    # check that extracted window is inside contig limit
    assert extract_contig_window(
        contig_size=16, positions_of_interest={15}, window_size=4
    ) == [(11, 15)]

    assert extract_contig_window(
        contig_size=10, positions_of_interest={2, 8}, window_size=2
    ) == [(0, 4), (6, 9)]

    # 12 window is (9,15)
    # 19  window is (16,22)
    # so when 12 and 19 are of interest window merge (9,22)
    assert extract_contig_window(
        contig_size=200, positions_of_interest={12}, window_size=3
    ) == [(9, 15)]
    assert extract_contig_window(
        contig_size=200, positions_of_interest={19}, window_size=3
    ) == [(16, 22)]
    assert extract_contig_window(
        contig_size=200, positions_of_interest={12, 19}, window_size=3
    ) == [(9, 22)]

    assert extract_contig_window(
        contig_size=10, positions_of_interest={2, 5, 8}, window_size=2
    ) == [(0, 9)]


def test_extract_contig_window_with_circular_contig():
    # TODO try to use @pytest.mark.parametrize to test different combinations
    # # check that circularity is properly taken into account
    assert extract_contig_window(
        contig_size=12, positions_of_interest={1}, window_size=2, is_circular=True
    ) == [(0, 3), (11, 11)]
    assert extract_contig_window(
        contig_size=12, positions_of_interest={1}, window_size=3, is_circular=True
    ) == [(0, 4), (10, 11)]
    assert extract_contig_window(
        contig_size=12, positions_of_interest={10}, window_size=3, is_circular=True
    ) == [(0, 1), (7, 11)]

    assert extract_contig_window(
        contig_size=12, positions_of_interest={6}, window_size=6, is_circular=True
    ) == [(0, 11)]
    assert extract_contig_window(
        contig_size=12, positions_of_interest={1}, window_size=6, is_circular=True
    ) == [(0, 11)]
    assert extract_contig_window(
        contig_size=12, positions_of_interest={1}, window_size=6, is_circular=False
    ) == [(0, 7)]

    assert extract_contig_window(
        contig_size=12, positions_of_interest={0, 9}, window_size=2, is_circular=False
    ) == [(0, 2), (7, 11)]

    assert extract_contig_window(
        contig_size=894,
        positions_of_interest=[151, 152, 153, 893],
        window_size=4,
        is_circular=True,
    ) == [(0, 3), (147, 157), (889, 893)]


def test_extract_contig_window_out_of_range():
    with pytest.raises(IndexError):
        extract_contig_window(contig_size=15, positions_of_interest={15}, window_size=1)

    with pytest.raises(IndexError):
        extract_contig_window(contig_size=15, positions_of_interest={-1}, window_size=1)


def test_get_n_next_genes_index():

    assert list(
        get_n_next_genes_index(
            current_index=6, next_genes_count=3, contig_size=100, is_circular=False
        )
    ) == [7, 8, 9]

    # there is no next gene because the current index is at the end of a non circular contig
    assert (
        list(
            get_n_next_genes_index(
                current_index=11, next_genes_count=2, contig_size=12, is_circular=False
            )
        )
        == []
    )


def test_get_n_next_genes_index_circular():
    assert list(
        get_n_next_genes_index(
            current_index=10, next_genes_count=3, contig_size=12, is_circular=True
        )
    ) == [11, 0, 1]
    assert list(
        get_n_next_genes_index(
            current_index=10, next_genes_count=16, contig_size=12, is_circular=True
        )
    ) == [11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]


def test_get_n_next_genes_index_out_of_range():
    with pytest.raises(IndexError):
        assert list(
            get_n_next_genes_index(
                current_index=10, next_genes_count=16, contig_size=8, is_circular=False
            )
        )


@pytest.fixture()
def simple_contig():

    contig = Contig(identifier=1, name="contig1", is_circular=False)

    contig_size = 6
    contig.length = contig_size
    genes = [Gene(str(i)) for i in range(contig_size)]
    organism = Organism("organism_A")
    for i, (gene, family_name) in enumerate(zip(genes, "ABCDEFGHIJKLMNOP")):
        family = GeneFamily(i, family_name)
        gene.fill_annotations(start=i + 1, stop=i + 2, strand="+", position=i)

        gene.fill_parents(organism, contig)

        contig.add(gene)
        family.add(gene)

    return contig


@pytest.fixture()
def simple_circular_contig():

    contig = Contig(identifier=2, name="contig2", is_circular=True)

    contig_size = 6
    genes = [Gene(str(i)) for i in range(contig_size)]

    for i, (gene, family_name) in enumerate(zip(genes, "ABCDEFGHIJKLMNOP")):
        family = GeneFamily(i, family_name)
        gene.fill_annotations(start=0, stop=0, strand=0, position=i)

        contig.add(gene)
        family.add(gene)

    return contig


def test_add_edges_to_context_graph(simple_contig):
    context_graph = nx.Graph()

    # simple_contig families : ABCDEF

    add_edges_to_context_graph(
        context_graph, contig=simple_contig, contig_windows=[(0, 3)], transitivity=1
    )

    nodes = sorted([n.name for n in context_graph.nodes()])
    edges = {tuple(sorted([n.name, v.name])) for n, v in context_graph.edges()}

    assert nodes == ["A", "B", "C", "D"]
    assert edges == {("A", "B"), ("A", "C"), ("B", "C"), ("B", "D"), ("C", "D")}


def test_add_edges_to_context_graph_2(simple_contig):
    context_graph = nx.Graph()

    # simple_contig families : A B-C-D E F

    add_edges_to_context_graph(
        context_graph, contig=simple_contig, contig_windows=[(1, 3)], transitivity=0
    )

    nodes = sorted([n.name for n in context_graph.nodes()])
    edges = {tuple(sorted([n.name, v.name])) for n, v in context_graph.edges()}

    assert nodes == ["B", "C", "D"]
    assert edges == {("B", "C"), ("C", "D")}


def test_add_edges_to_context_graph_linear(simple_contig):

    #    genes : 1-2-3-4-5-6
    # families : A-B-C-D-E-F
    #  windows : _____   ___ [(0,2) (4,5)]

    context_graph = nx.Graph()

    add_edges_to_context_graph(
        context_graph,
        contig=simple_contig,
        contig_windows=[(4, 5), (0, 2)],
        transitivity=0,
    )

    nodes = sorted([n.name for n in context_graph.nodes()])
    edges = {tuple(sorted([n.name, v.name])) for n, v in context_graph.edges()}

    assert nodes == ["A", "B", "C", "E", "F"]
    assert edges == {
        ("A", "B"),
        ("B", "C"),
        ("E", "F"),
    }


def test_add_edges_to_context_graph_circular(simple_contig):

    #    genes : 1-2-3-4-5-6
    # families : A-B-C-D-E-F
    #  windows : _____   ___ [(0,2) (4,5)]

    context_graph = nx.Graph()
    simple_contig.is_circular = True
    add_edges_to_context_graph(
        context_graph,
        contig=simple_contig,
        contig_windows=[(4, 5), (0, 2)],
        transitivity=0,
    )

    nodes = sorted([n.name for n in context_graph.nodes()])
    edges = {tuple(sorted([n.name, v.name])) for n, v in context_graph.edges()}

    assert nodes == ["A", "B", "C", "E", "F"]
    assert edges == {
        ("A", "B"),
        ("B", "C"),
        ("E", "F"),
        ("A", "F"),
    }  # circular so F and A are linked


def test_compute_gene_context_graph(simple_contig):

    #              genes : 0-1-2-3-4-5
    #           families : A-B-C-D-E-F
    # family of interest :     ^
    #       windows of 2 : ___   ___

    # simple case with only one contig with 6 genes and 6 families

    families_in_contigs = [g.family for g in simple_contig.genes]
    family_names_of_interest = ["C"]
    families_of_interest = {
        f for f in families_in_contigs if f.name in family_names_of_interest
    }

    context_graph, _ = compute_gene_context_graph(
        families_of_interest, transitive=0, window_size=2
    )
    nodes = sorted([n.name for n in context_graph.nodes()])
    edges = {tuple(sorted([n.name, v.name])) for n, v in context_graph.edges()}

    assert nodes == ["A", "B", "C", "D", "E"]
    assert edges == {("A", "B"), ("B", "C"), ("C", "D"), ("D", "E")}
