#!/usr/bin/env python3
#coding: utf8

#default libraries
from collections import defaultdict


class Edge:
    """The Edge class represents an edge between two gene families in the pangenome graph. It is associated with all the organisms in which the neighborship is found, and all the involved genes as well.

    :param sourceGene: a first gene to initialize the edge
    :type sourceGene: :class:`ppanggolin.genome.Gene`
    :param targetGene: a second gene to initialize the edge
    :type targetGene: :class:`ppanggolin.genome.Gene`
    """
    def __init__(self, sourceGene, targetGene):
        if sourceGene.family is None:
            raise Exception(f"You cannot create a graph without gene families. gene {sourceGene.ID} did not have a gene family.")
        if targetGene.family is None:
            raise Exception(f"You cannot create a graph without gene families. gene {targetGene.ID} did not have a gene family.")
        self.source = sourceGene.family
        self.target = targetGene.family
        self.source._edges[self.target] = self
        self.target._edges[self.source] = self
        self.organisms = defaultdict(list)
        self.addGenes(sourceGene, targetGene)

    def getOrgDict(self):
        """
        
        :return: A dictionnary of the Organisms in which the edge is found, with organisms as key and an iterable of the pairs of genes as value
        :rtype: dict[:class:`ppanggolin.genome.Organism`, tuple[:class:`ppanggolin.genome.Gene`, :class:`ppanggolin.genome.Gene`]]
        """
        return self.organisms

    @property
    def genePairs(self):
        """
        
        :return: A list of all the gene pairs of the Edge
        :rtype: list[tuple[:class:`ppanggolin.genome.Gene`, :class:`ppanggolin.genome.Gene`]]
        """
        return [ gene_pair for gene_list in self.organisms.values() for gene_pair in gene_list ]

    def addGenes(self, sourceGene, targetGene):
        """Adds genes to the edge. They are supposed to be on the same organism.

        :param sourceGene: a source gene to add to the edge
        :type sourceGene: :class:`ppanggolin.genome.Gene`
        :param targetGene: a target gene to add to the edge
        :type targetGene: :class:`ppanggolin.genome.Gene`
        :raises Exception: If the genes are not on the same organism.
        """
        org = sourceGene.organism
        if org != targetGene.organism:
            raise Exception(f"You tried to create an edge between two genes that are not even in the same organism ! (genes are '{sourceGene.ID}' and '{targetGene.ID}')")
        self.organisms[org].append((sourceGene, targetGene))