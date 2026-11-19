#! /usr/bin/env python3

from typing import Dict, List, Optional, Set, Tuple

import pytest

from ppanggolin.genome import Contig, Gene, Organism
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.graph.makeGraph import add_contig_edges
from ppanggolin.pangenome import Pangenome


def build_contig(
    families: List[str],
    fragments: Optional[Set[int]] = None,
    removed: Optional[Set[str]] = None,
    is_circular: bool = False,
) -> Tuple[Pangenome, Contig]:
    """Build a pangenome with one contig whose genes belong to the given families.

    :param families: Family name of each gene, in contig order
    :param fragments: Positions of the genes flagged as fragments
    :param removed: Names of the families removed from the graph
    :param is_circular: Whether the contig is circular
    :return: The pangenome and its contig
    """
    fragments = fragments or set()
    removed = removed or set()
    pangenome = Pangenome()
    organism = Organism("organism")
    contig = Contig(0, "contig", is_circular=is_circular)
    organism.add(contig)
    fams: Dict[str, GeneFamily] = {}
    for position, fam_name in enumerate(families):
        if fam_name not in fams:
            fams[fam_name] = GeneFamily(len(fams), fam_name)
            fams[fam_name].removed = fam_name in removed
            pangenome.add_gene_family(fams[fam_name])
        gene = Gene(f"gene_{position}")
        gene.fill_annotations(
            start=10 * position + 1,
            stop=10 * position + 9,
            strand="+",
            position=position,
        )
        gene.is_fragment = position in fragments
        contig.add(gene)
        gene.fill_parents(organism, contig)
        fams[fam_name].add(gene)
    return pangenome, contig


def gene_pairs(pangenome: Pangenome) -> Set[frozenset]:
    """Unordered pairs of gene identifiers linked in the graph."""
    return {
        frozenset((gene1.ID, gene2.ID))
        for edge in pangenome.edges
        for gene1, gene2 in edge.gene_pairs
    }


def pair(position1: int, position2: int) -> frozenset:
    return frozenset((f"gene_{position1}", f"gene_{position2}"))


@pytest.mark.parametrize("is_circular", [False, True])
def test_consecutive_genes_are_linked(is_circular):
    pangenome, contig = build_contig(["A", "B", "C"], is_circular=is_circular)
    add_contig_edges(pangenome, contig)
    expected = {pair(0, 1), pair(1, 2)} | ({pair(2, 0)} if is_circular else set())
    assert gene_pairs(pangenome) == expected


def test_fragments_of_a_same_family_are_not_linked():
    pangenome, contig = build_contig(["A", "B", "B", "C"], fragments={2})
    add_contig_edges(pangenome, contig)
    assert gene_pairs(pangenome) == {pair(0, 1), pair(2, 3)}


def test_tandem_copies_without_fragment_are_linked():
    pangenome, contig = build_contig(["A", "B", "B", "C"])
    add_contig_edges(pangenome, contig)
    assert pair(1, 2) in gene_pairs(pangenome)


def test_circular_closure_skips_fragments_of_a_same_family():
    # first and last genes are pieces of the same fragmented gene, split by the contig origin
    pangenome, contig = build_contig(
        ["A", "B", "C", "A"], fragments={0, 3}, is_circular=True
    )
    add_contig_edges(pangenome, contig)
    assert gene_pairs(pangenome) == {pair(0, 1), pair(1, 2), pair(2, 3)}


def test_circular_closure_skips_removed_families():
    pangenome, contig = build_contig(
        ["R", "A", "B", "C"], removed={"R"}, is_circular=True
    )
    add_contig_edges(pangenome, contig)
    assert gene_pairs(pangenome) == {pair(1, 2), pair(2, 3), pair(3, 1)}


def test_removed_families_are_skipped():
    pangenome, contig = build_contig(["A", "R", "B"], removed={"R"})
    add_contig_edges(pangenome, contig)
    assert gene_pairs(pangenome) == {pair(0, 2)}


def test_circular_contig_with_a_single_gene_has_no_edge():
    pangenome, contig = build_contig(["A"], is_circular=True)
    add_contig_edges(pangenome, contig)
    assert gene_pairs(pangenome) == set()


def test_circular_contig_with_two_genes_is_linked_twice():
    pangenome, contig = build_contig(["A", "B"], is_circular=True)
    add_contig_edges(pangenome, contig)
    (edge,) = pangenome.edges
    assert len(edge.gene_pairs) == 2
