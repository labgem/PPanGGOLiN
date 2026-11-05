#!/usr/bin/env python3

# default libraries
from collections.abc import Iterable
from dataclasses import dataclass

# installed libraries
import networkx as nx
import numpy as np
import scipy.sparse as sp

# local libraries
from ppanggolin.genome import Organism
from ppanggolin.pangenome import Pangenome


@dataclass
class NemIndex:
    """NEM's input for a set of genomes, before it is made dense.

    :param org_names: Genome names, one per column.
    :param fam_names: Gene family names, one per row of `presence`.
    :param presence: Families x genomes presence matrix.
    :param coverage: Edges x genomes gene-pair counts.
    :param edge_src: Row of each edge's source family.
    :param edge_tgt: Row of each edge's target family.
    """

    org_names: list[str]
    fam_names: list[str]
    presence: sp.csr_matrix
    coverage: sp.csr_matrix
    edge_src: np.ndarray
    edge_tgt: np.ndarray


@dataclass
class NemInput:
    """Dense NEM input for one set of genomes.

    :param presence: Families x genomes presence matrix.
    :param graph: Smoothing graph between family rows, edges carrying `weight`.
    :param fam_names: Gene family names, one per row of `presence`.
    :param edges_weight: Total weight of the smoothing edges.
    """

    presence: np.ndarray
    graph: nx.DiGraph
    fam_names: list[str]
    edges_weight: float

    @property
    def nb_families(self) -> int:
        """Number of gene families."""
        return self.presence.shape[0]

    @property
    def nb_genomes(self) -> int:
        """Number of genomes."""
        return self.presence.shape[1]

    def scaled_beta(self, beta: float) -> float:
        """Scale the smoothing strength to the size of the graph.

        :param beta: Smoothing strength given by the user.
        :return: The strength NEM is run with.
        """
        return beta * (self.nb_families / self.edges_weight)


def ordered_organisms(organisms: Iterable[Organism]) -> list[Organism]:
    """Sort genomes in a stable, process-independent order.
    :param organisms: Genomes to order.
    :return: The genomes sorted by name.
    """
    return sorted(organisms, key=lambda org: org.name)


def _presence_matrix(families: list, org_index: dict) -> sp.csr_matrix:
    """Build the families x genomes presence matrix.

    :param families: Gene families, one per row.
    :param org_index: Column of each genome.
    :return: The presence matrix.
    """
    nb_presences = sum(fam.number_of_organisms for fam in families)
    rows = np.empty(nb_presences, dtype=np.int32)
    cols = np.empty(nb_presences, dtype=np.int32)
    pos = 0
    for i, fam in enumerate(families):
        for org in fam.organisms:
            rows[pos] = i
            cols[pos] = org_index[org]
            pos += 1
    return sp.csr_matrix(
        (np.ones(nb_presences, dtype=np.int8), (rows, cols)),
        shape=(len(families), len(org_index)),
    )


def _coverage_matrix(
    pangenome: Pangenome, fam_index: dict, org_index: dict
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray]:
    """Count the gene pairs of each edge in each genome.

    :param pangenome: Pangenome with its neighbors graph.
    :param fam_index: Row of each gene family.
    :param org_index: Column of each genome.
    :return: The edges x genomes counts, and the source and target row of each edge.
    """
    nb_edges = pangenome.number_of_edges
    nb_covered = sum(edge.number_of_organisms for edge in pangenome.edges)
    src = np.empty(nb_edges, dtype=np.int32)
    tgt = np.empty(nb_edges, dtype=np.int32)
    rows = np.empty(nb_covered, dtype=np.int32)
    cols = np.empty(nb_covered, dtype=np.int32)
    counts = np.empty(nb_covered, dtype=np.int32)
    pos = 0
    for e, edge in enumerate(pangenome.edges):
        src[e] = fam_index[edge.source]
        tgt[e] = fam_index[edge.target]
        for org, gene_pairs in edge.get_organisms_dict().items():
            rows[pos] = e
            cols[pos] = org_index[org]
            counts[pos] = len(gene_pairs)
            pos += 1
    coverage = sp.csr_matrix((counts, (rows, cols)), shape=(nb_edges, len(org_index)))
    return coverage, src, tgt


def build_nem_index(pangenome: Pangenome) -> NemIndex:
    """Precompute NEM's input for all the genomes of a pangenome.

    :param pangenome: Pangenome with gene families and neighbors graph.
    :return: The index.
    """
    organisms = ordered_organisms(pangenome.organisms)
    org_index = {org: i for i, org in enumerate(organisms)}
    families = list(pangenome.gene_families)
    fam_index = {fam: i for i, fam in enumerate(families)}

    presence = _presence_matrix(families, org_index)
    coverage, src, tgt = _coverage_matrix(pangenome, fam_index, org_index)

    return NemIndex(
        org_names=[org.name for org in organisms],
        fam_names=[fam.name for fam in families],
        presence=presence,
        coverage=coverage,
        edge_src=src,
        edge_tgt=tgt,
    )


def subset_nem_index(index: NemIndex, genomes: Iterable[str]) -> NemIndex:
    """Restrict an index to some of its genomes.

    :param index: Index to slice.
    :param genomes: Names of the genomes to keep.
    :return: The index of these genomes, without the families and edges they do not have.
    """
    genomes = list(genomes)
    column = {name: i for i, name in enumerate(index.org_names)}
    cols = np.sort(
        np.fromiter(
            (column[name] for name in genomes), dtype=np.int32, count=len(genomes)
        )
    )

    presence = index.presence[:, cols]
    rows_kept = np.flatnonzero(np.diff(presence.indptr))
    local_row = np.full(presence.shape[0], -1, dtype=np.int32)
    local_row[rows_kept] = np.arange(rows_kept.size, dtype=np.int32)

    coverage = index.coverage[:, cols]
    live = np.flatnonzero(np.diff(coverage.indptr))

    return NemIndex(
        org_names=[index.org_names[c] for c in cols],
        fam_names=[index.fam_names[r] for r in rows_kept],
        presence=presence[rows_kept],
        coverage=coverage[live],
        edge_src=local_row[index.edge_src[live]],
        edge_tgt=local_row[index.edge_tgt[live]],
    )


def _smoothing_graph(
    index: NemIndex, max_degree_smoothing: int
) -> tuple[nx.DiGraph, float]:
    """Link each family to its neighbors, for the families NEM smooths.

    :param index: Index of the genomes to partition.
    :param max_degree_smoothing: Maximum degree of a family included in the smoothing.
    :return: The graph and the total weight of its edges.
    """
    nb_families, nb_genomes = index.presence.shape
    src, tgt = index.edge_src, index.edge_tgt
    score = np.asarray(index.coverage.sum(axis=1)).ravel() / nb_genomes

    other = src != tgt
    degree = np.bincount(src, minlength=nb_families) + np.bincount(
        tgt[other], minlength=nb_families
    )
    sum_dist = np.bincount(src, weights=score, minlength=nb_families) + np.bincount(
        tgt[other], weights=score[other], minlength=nb_families
    )

    smoothed = (degree > 0) & (degree < max_degree_smoothing)
    total_edges_weight = sum_dist[smoothed].sum()

    from_src, from_tgt = smoothed[src], smoothed[tgt] & other
    graph = nx.DiGraph()
    graph.add_nodes_from(range(nb_families))
    graph.add_weighted_edges_from(
        zip(
            np.concatenate([src[from_src], tgt[from_tgt]]).tolist(),
            np.concatenate([tgt[from_src], src[from_tgt]]).tolist(),
            np.concatenate([score[from_src], score[from_tgt]]).tolist(),
        )
    )
    return graph, total_edges_weight / 2


def build_nem_input(index: NemIndex, max_degree_smoothing: int) -> NemInput:
    """Build NEM's input for all the genomes of an index.

    :param index: Index of the genomes to partition.
    :param max_degree_smoothing: Maximum degree of a family included in the smoothing.
    :return: The dense input.
    """
    graph, edges_weight = _smoothing_graph(index, max_degree_smoothing)
    return NemInput(
        presence=index.presence.toarray().astype(np.float64),
        graph=graph,
        fam_names=index.fam_names,
        edges_weight=edges_weight,
    )
