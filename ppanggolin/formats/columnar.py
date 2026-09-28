from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
import tables


def locate(gene_ids: np.ndarray, columns: Sequence[np.ndarray]) -> List[np.ndarray]:
    """
    For each column of gene identifiers, the row it occupies in ``gene_ids``.

    :param gene_ids: gene identifiers, in the order to resolve against
    :param columns: arrays of gene identifiers to locate

    :return: one array of row indices per column
    """
    sizes = [gene_ids.size] + [column.size for column in columns]
    bounds = np.cumsum(sizes)
    _, inverse = np.unique(np.concatenate([gene_ids, *columns]), return_inverse=True)
    slot = np.full(int(inverse.max()) + 1 if inverse.size else 0, -1, dtype=np.int64)
    slot[inverse[: gene_ids.size]] = np.arange(gene_ids.size)

    located = []
    for start, stop in zip(bounds[:-1], bounds[1:]):
        rows = slot[inverse[start:stop]]
        if rows.size and rows.min() < 0:
            raise KeyError(
                "The pangenome references genes absent from the annotations."
            )
        located.append(rows)
    return located


def gene_genome_arrays(h5f: tables.File) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Gene identifiers, the genome index of each gene, and the genome names.

    :param h5f: open pangenome file

    :return: gene identifiers, genome index per gene, genome names
    """
    contigs = h5f.root.annotations.contigs
    contig_ids = contigs.read(field="ID")
    genomes, contig_genome = np.unique(
        contigs.read(field="genome"), return_inverse=True
    )

    genes = h5f.root.annotations.genes
    gene_ids = genes.read(field="ID")
    if contig_ids.size == 0:
        return gene_ids, np.empty(0, dtype=np.int64), genomes

    # Contig ids are small dense integers, so a lookup array beats a dict.
    lookup = np.empty(contig_ids.max() + 1, dtype=np.int64)
    lookup[contig_ids] = contig_genome

    return gene_ids, lookup[genes.read(field="contig")], genomes


def _by_first_appearance(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Factorise ``values``, numbering the codes by first appearance.

    :param values: column to factorise

    :return: the distinct values in first-appearance order, and a code per row
    """
    uniques, first, inverse = np.unique(values, return_index=True, return_inverse=True)
    order = np.argsort(first)
    rank = np.empty(order.size, dtype=np.int64)
    rank[order] = np.arange(order.size)
    return uniques[order], rank[inverse]


def _pair_counts(left: np.ndarray, right: np.ndarray, width: int):
    """Distinct ``(left, right)`` pairs and how many rows each covers."""
    keys, counts = np.unique(left * width + right, return_counts=True)
    rows, cols = np.divmod(keys, width)
    return rows, cols, counts


def build_index(h5f: tables.File) -> Dict:
    """
    Everything ``partition`` needs from a pangenome, without the object graph.

    :param h5f: open pangenome file, clustered and with a neighbours graph

    :return: family names, genome names, the families x genomes presence matrix,
             the edges x genomes coverage matrix and the edge endpoints
    """
    gene_ids, genome_of_gene, genomes = gene_genome_arrays(h5f)
    n_genomes = genomes.size

    families = h5f.root.geneFamilies
    fam_names, fam_of_row = _by_first_appearance(families.read(field="geneFam"))
    edges = h5f.root.edges

    family_genes, source_genes, target_genes = locate(
        gene_ids,
        [
            families.read(field="gene"),
            edges.read(field="geneSource"),
            edges.read(field="geneTarget"),
        ],
    )

    # A family is present in a genome if any of its genes is.
    rows, cols, _ = _pair_counts(fam_of_row, genome_of_gene[family_genes], n_genomes)
    presence = sp.csr_matrix(
        (np.ones(rows.size, dtype=np.int8), (rows, cols)),
        shape=(fam_names.size, n_genomes),
    )
    del rows, cols

    # An edge joins two families whichever way round the adjacency was read, so
    # it is keyed on the unordered pair; its orientation is that of the first
    # adjacency that produced it, which is what `Pangenome.add_edge` records.
    row_to_family = np.empty(gene_ids.size, dtype=np.int64)
    row_to_family[family_genes] = _by_first_appearance(families.read(field="geneFam"))[
        1
    ]
    source_family = row_to_family[source_genes]
    target_family = row_to_family[target_genes]
    low = np.minimum(source_family, target_family)
    high = np.maximum(source_family, target_family)

    _, first, edge_of_row = np.unique(
        low * fam_names.size + high, return_index=True, return_inverse=True
    )
    order = np.argsort(first)
    rank = np.empty(order.size, dtype=np.int64)
    rank[order] = np.arange(order.size)
    edge_of_row = rank[edge_of_row]
    first = first[order]

    edge_rows, edge_cols, counts = _pair_counts(
        edge_of_row, genome_of_gene[source_genes], n_genomes
    )
    coverage = sp.csr_matrix(
        (counts.astype(np.int32), (edge_rows, edge_cols)),
        shape=(first.size, n_genomes),
    )

    return {
        "org_index": {name.decode(): i for i, name in enumerate(genomes)},
        "fam_names": [name.decode() for name in fam_names],
        "presence": presence,
        "coverage": coverage,
        "edge_src": source_family[first].astype(np.int32),
        "edge_tgt": target_family[first].astype(np.int32),
    }
