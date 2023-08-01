#!/usr/bin/env python3
# coding: utf8

# default libraries
from collections import defaultdict
from typing import Dict, List, Tuple

from ppanggolin.genome import Gene, Organism


class Edge:
    """The Edge class represents an edge between two gene families in the pangenome graph. It is associated with all the
     organisms in which the neighborship is found, and all the involved genes as well.

    :param source_gene: a first gene to initialize the edge
    :param target_gene: a second gene to initialize the edge
    """

    def __init__(self, source_gene: Gene, target_gene: Gene):
        if source_gene.family is None:
            raise Exception(f"You cannot create a graph without gene families. "
                            f"gene {source_gene.ID} did not have a gene family.")
        if target_gene.family is None:
            raise Exception(f"You cannot create a graph without gene families. "
                            f"gene {target_gene.ID} did not have a gene family.")
        self.source = source_gene.family
        self.target = target_gene.family
        self.source.set_edge(self.target, self)
        self.target.set_edge(self.source, self)
        self.organisms = defaultdict(list)
        self.add_genes(source_gene, target_gene)

    def get_org_dict(self) -> Dict[Organism, List[Tuple[Gene, Gene]]]:
        """ Create a dictionnary of the Organisms in which the edge is found

        :return: Dictionary with organisms as key and an iterable of the pairs of genes as value
        """
        return self.organisms

    @property
    def gene_pairs(self) -> List[Tuple[Gene, Gene]]:
        """ Get list of all the gene pairs of the Edge

        :return: A list of all the gene pairs of the Edge
        """
        return [gene_pair for gene_list in self.organisms.values() for gene_pair in gene_list]

    def add_genes(self, source_gene: Gene, target_gene: Gene):
        """Adds genes to the edge. They are supposed to be on the same organism.

        :param source_gene: a source gene to add to the edge
        :param target_gene: a target gene to add to the edge

        :raises Exception: If the genes are not on the same organism.
        """
        org = source_gene.organism
        if org != target_gene.organism:
            raise Exception(f"You tried to create an edge between two genes that are not even in the same organism ! "
                            f"(genes are '{source_gene.ID}' and '{target_gene.ID}')")
        self.organisms[org].append((source_gene, target_gene))
