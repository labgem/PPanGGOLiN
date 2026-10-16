#!/usr/bin/env python3

# default libraries
from dataclasses import dataclass

# installed libraries
import networkx as nx
import numpy as np
import scipy.sparse as sp

# local libraries
from ppanggolin.pangenome import Pangenome


@dataclass
class NemIndex:
    """
    NEM's input for a set of genomes, before it is made dense

    :param org_names: genome names, one per column
    :param fam_names: gene family names, one per row of `presence`
    :param presence: families x genomes presence matrix
    :param coverage: edges x genomes gene-pair counts
    :param edge_src: row of each edge's source family
    :param edge_tgt: row of each edge's target family
    """

    org_names: list
    fam_names: list
    presence: sp.csr_matrix
    coverage: sp.csr_matrix
    edge_src: np.ndarray
    edge_tgt: np.ndarray


def ordered_organisms(organisms) -> list:
    """
    Genomes in a stable, process-independent order.
    :param organisms: genomes to order
    :return: the genomes sorted by name
    """
    return sorted(organisms, key=lambda org: org.name)


def build_nem_index(pangenome: Pangenome) -> NemIndex:
    """
    Precompute NEM's input for all the genomes of a pangenome

    :param pangenome: Pangenome with gene families and graph

    :return: the index
    """
    organisms = ordered_organisms(pangenome.organisms)
    org_index = {org: i for i, org in enumerate(organisms)}
    families = list(pangenome.gene_families)
    fam_index = {fam: i for i, fam in enumerate(families)}

    nb_presences = sum(fam.number_of_organisms for fam in families)
    rows = np.empty(nb_presences, dtype=np.int32)
    cols = np.empty(nb_presences, dtype=np.int32)
    pos = 0
    for i, fam in enumerate(families):
        for org in fam.organisms:
            rows[pos] = i
            cols[pos] = org_index[org]
            pos += 1
    presence = sp.csr_matrix(
        (np.ones(nb_presences, dtype=np.int8), (rows, cols)),
        shape=(len(families), len(organisms)),
    )
    del rows, cols

    nb_edges = pangenome.number_of_edges
    nb_covered = sum(edge.number_of_organisms for edge in pangenome.edges)
    src = np.empty(nb_edges, dtype=np.int32)
    tgt = np.empty(nb_edges, dtype=np.int32)
    e_rows = np.empty(nb_covered, dtype=np.int32)
    e_cols = np.empty(nb_covered, dtype=np.int32)
    counts = np.empty(nb_covered, dtype=np.int32)
    pos = 0
    for e, edge in enumerate(pangenome.edges):
        src[e] = fam_index[edge.source]
        tgt[e] = fam_index[edge.target]
        for org, gene_pairs in edge.get_organisms_dict().items():
            e_rows[pos] = e
            e_cols[pos] = org_index[org]
            counts[pos] = len(gene_pairs)
            pos += 1
    coverage = sp.csr_matrix(
        (counts, (e_rows, e_cols)), shape=(nb_edges, len(organisms))
    )
    del e_rows, e_cols, counts

    return NemIndex(
        org_names=[org.name for org in organisms],
        fam_names=[fam.name for fam in families],
        presence=presence,
        coverage=coverage,
        edge_src=src,
        edge_tgt=tgt,
    )


def subset_nem_index(index: NemIndex, genomes) -> NemIndex:
    """
    Restrict an index to some of its genomes

    :param index: index to slice
    :param genomes: names of the genomes to keep

    :return: the index of these genomes, without the families and edges they do not have
    """
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


def build_nem_input(index: NemIndex, max_degree_smoothing: int = 10) -> tuple:
    """
    Build NEM's inputs for all the genomes of an index

    :param index: index of the genomes to partition
    :param max_degree_smoothing: maximum degree of a node included in the smoothing

    :return: (presence, graph, fam_names, edges_weight, nb_families)
    """
    nb_families, nb_genomes = index.presence.shape
    presence = index.presence.toarray().astype(np.float64)

    src, tgt = index.edge_src, index.edge_tgt
    score = np.asarray(index.coverage.sum(axis=1)).ravel() / nb_genomes

    loop = src == tgt
    other = ~loop
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

    return presence, graph, index.fam_names, total_edges_weight / 2, nb_families
