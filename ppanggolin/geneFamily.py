#!/usr/bin/env python3
#coding: utf8

#default libraries
from collections import defaultdict

#installed libraries
import gmpy2

#local libraries
from ppanggolin.genome import Gene


class GeneFamily:
    def __init__(self, ID, name):
        self.name = name
        self.ID = ID
        self._edges = {}
        self._genePerOrg = defaultdict(set)
        self.genes = set()
        self.removed = False#for the repeated family not added in the main graph
        self.sequence = ""
        self.partition = ""

    def addSequence(self, seq):
        self.sequence = seq

    def addPartition(self, partition):
        self.partition = partition

    @property
    def namedPartition(self):
        if self.partition == "":
            raise Exception("The gene family has not beed associated to a partition")
        if self.partition.startswith("P"):
            return "persistent"
        elif self.partition.startswith("C"):
            return "cloud"
        elif self.partition.startswith("S"):
            return "shell"
        else:
            return "undefined"

    def addGene(self, gene):
        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type object was expected, but '{type(gene)}' type object was provided.")
        self.genes.add(gene)
        gene.family = self
        if hasattr(gene, "organism"):
            self._genePerOrg[gene.organism].add(gene)

    def mkBitarray(self, index):
        """ produces a bitarray representing the presence / absence of the family in the pangenome"""
        self.bitarray = gmpy2.xmpz(0)
        for org in self.organisms:
            self.bitarray[index[org]] = 1

    def getOrgDict(self):
        try:
            return self._genePerOrg
        except AttributeError:
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg

    def getGenesPerOrg(self, org):
        try:
            return self._genePerOrg[org]
        except AttributeError:
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg[org]

    @property
    def neighbors(self):
        return set(self._edges.keys())

    @property
    def edges(self):
        return self._edges.values()

    @property
    def organisms(self):
        try:
            return self._genePerOrg.keys()
        except AttributeError:#then the genes have been added before they had organisms
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg.keys()
