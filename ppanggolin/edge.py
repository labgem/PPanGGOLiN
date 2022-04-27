#!/usr/bin/env python3
# coding: utf8

# default libraries
from collections import defaultdict


class Edge:
    """The Edge class represents an edge between two gene families in the pan graph. It is associated with all the
     organisms in which the neighborship is found, and all the involved genes as well.

    :param source_gene: a first gene to initialize the edge
    :type source_gene: :class:`ppanggolin.genome.Gene`
    :param target_gene: a second gene to initialize the edge
    :type target_gene: :class:`ppanggolin.genome.Gene`
    """

    def __init__(self, source_gene, target_gene):
        if source_gene.family is None:
            raise Exception(
                f"You cannot create a graph without gene families. gene {source_gene.ID} did not have a gene family.")
        if target_gene.family is None:
            raise Exception(
                f"You cannot create a graph without gene families. gene {target_gene.ID} did not have a gene family.")
        self.source = source_gene.family
        self.target = target_gene.family
        self.source._edges[self.target] = self
        self.target._edges[self.source] = self
        self.organisms = defaultdict(list)
        self.add_genes(source_gene, target_gene)

    def get_org_dict(self):
        """

        :return: A dictionnary of the Organisms in which the edge is found, with organisms as key and an iterable of the
         pairs of genes as value
        :rtype: dict[:class:`ppanggolin.genome.Organism`, tuple[:class:`ppanggolin.genome.Gene`,
        :class:`ppanggolin.genome.Gene`]]
        """
        return self.organisms

    @property
    def gene_pairs(self):
        """
        
        :return: A list of all the gene pairs of the Edge
        :rtype: list[tuple[:class:`ppanggolin.genome.Gene`, :class:`ppanggolin.genome.Gene`]]
        """
        return [gene_pair for gene_list in self.organisms.values() for gene_pair in gene_list]

    def add_genes(self, source_gene, target_gene):
        """Adds genes to the edge. They are supposed to be on the same organism.

        :param source_gene: a source gene to add to the edge
        :type source_gene: :class:`ppanggolin.genome.Gene`
        :param target_gene: a target gene to add to the edge
        :type target_gene: :class:`ppanggolin.genome.Gene`
        :raises Exception: If the genes are not on the same organism.
        """
        org = source_gene.organism
        if org != target_gene.organism:
            raise Exception(f"You tried to create an edge between two genes that are not even in the same organism ! "
                            f"(genes are '{source_gene.ID}' and '{target_gene.ID}')")
        self.organisms[org].append((source_gene, target_gene))
