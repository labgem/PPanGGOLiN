#!/usr/bin/env python3
#coding: utf8

#default libraries
from collections import defaultdict

#installed libraries
import gmpy2

#local libraries
from ppanggolin.genome import Gene


class GeneFamily:
    """This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.

    """
    def __init__(self, ID, name):
        """Constructor method

        :param ID: The internal identifier to give to the gene family
        :type ID: any
        :param name: The name of the gene family (to be printed in output files)
        :type name: str
        """
        self.name = str(name)
        self.ID = ID
        self._edges = {}
        self._genePerOrg = defaultdict(set)
        self.genes = set()
        self.removed = False#for the repeated family not added in the main graph
        self.sequence = ""
        self.partition = ""

    def addSequence(self, seq):
        """Assigns a protein sequence to the gene family.

        :param seq: the sequence to add to the gene family
        :type seq: str
        """
        self.sequence = seq

    def addPartition(self, partition):
        """Assigns a partition to the gene family. It should be the raw partition name provided by NEM.

        :param partition: The partition
        :type partition: str
        """
        self.partition = partition

    @property
    def namedPartition(self):
        """Reads the :attr:partition attribute and returns a meaningful name

        :raises Exception: If the gene family has no partition assigned
        :return: the partition name of the gene family
        :rtype: str
        """
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
        """Add a gene to the gene family, and sets the gene's :attr:family accordingly.

        :param gene: the gene to add
        :type gene: :class:`ppanggolin.genome.Gene`
        :raises TypeError: If the provided `gene` is of the wrong type
        """
        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type object was expected, but '{type(gene)}' type object was provided.")
        self.genes.add(gene)
        gene.family = self
        if hasattr(gene, "organism"):
            self._genePerOrg[gene.organism].add(gene)

    def mkBitarray(self, index):
        """Produces a bitarray representing the presence / absence of the family in the pangenome using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.

        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`
        :type index: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        self.bitarray = gmpy2.xmpz(0)#pylint: disable=no-member
        for org in self.organisms:
            self.bitarray[index[org]] = 1

    def getOrgDict(self):
        """Returns the organisms and the genes belonging to the gene family

        :return: a dictionnary of organism as key and set of genes as values
        :rtype: dict[ :class:`ppanggolin.genome.Organism` ,set[:class:`ppanggolin.genome.Gene`]
        """
        try:
            return self._genePerOrg
        except AttributeError:
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg

    def getGenesPerOrg(self, org):
        """Returns the genes belonging to the gene family in the given Organism

        :param org: Organism to look for
        :type org: :class:`ppanggolin.genome.Organism`
        :return: a set of gene(s)
        :rtype: set[:class:`ppanggolin.genome.Gene`]
        """
        try:
            return self._genePerOrg[org]
        except AttributeError:
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg[org]

    @property
    def neighbors(self):
        """Returns all of the :class:`ppanggolin.geneFamily.GeneFamily` that are linked with an edge

        :return: Neighbors
        :rtype: set[:class:`ppanggolin.geneFamily.GeneFamily`]
        """
        return set(self._edges.keys())

    @property
    def edges(self):
        """Returns all of the :class:`ppanggolin.pangenome.Edge` that are linked to this gene family

        :return: Edges of the gene family
        :rtype: list[:class:`ppanggolin.pangenome.Edge`]
        """
        return list(self._edges.values())

    @property
    def organisms(self):
        """Returns all of the :class:`ppanggolin.genome.Organism` that have this gene family

        :return: Organisms that have this gene family
        :rtype: set[:class:`ppanggolin.genome.Organism`]
        """
        try:
            return set(self._genePerOrg.keys())
        except AttributeError:#then the genes have been added before they had organisms
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return set(self._genePerOrg.keys())