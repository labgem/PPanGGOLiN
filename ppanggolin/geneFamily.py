#!/usr/bin/env python3
# coding: utf8

# default libraries
from __future__ import annotations
from collections import defaultdict
import logging

# installed libraries
from typing import Dict, List, Set

import gmpy2

# local libraries
from ppanggolin.edge import Edge
from ppanggolin.genome import Gene, Organism
from ppanggolin.metadata import MetaFeatures


class GeneFamily(MetaFeatures):
    """
    This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.

    :param family_id: The internal identifier to give to the gene family
    :type family_id: any
    :param name: The name of the gene family (to be printed in output files)
    :type name: str
    """

    def __init__(self, family_id: int, name: str):
        assert isinstance(family_id, int), "GeneFamily object id should be an integer"
        assert isinstance(name, str), "GeneFamily object name should be a string"
        assert name != '', "GeneFamily object cannot be created with an empty name"

        super().__init__()
        self.name = str(name)
        self.ID = family_id
        self._edges = {}
        self._genePerOrg = defaultdict(set)
        self.genes = set()
        self.removed = False  # for the repeated family not added in the main graph
        self.sequence = ""
        self.partition = ""
        self.spot = set()
        self.modules = set()
        self.bitarray = None

    def add_sequence(self, seq: str):
        """Assigns a protein sequence to the gene family.

        :param seq: the sequence to add to the gene family
        """
        assert isinstance(seq, str) and str != "", "Sequence must be a string and not empty"
        self.sequence = seq

    def add_partition(self, partition: str):
        """Assigns a partition to the gene family. It should be the raw partition name provided by NEM.

        :param partition: The partition
        """
        self.partition = partition

    @property
    def named_partition(self) -> str:
        """Reads the partition attribute and returns a meaningful name

        :raises Exception: If the gene family has no partition assigned

        :return: the partition name of the gene family
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

    def add_gene(self, gene: Gene):
        """Add a gene to the gene family, and sets the gene's :attr:family accordingly.

        :param gene: the gene to add

        :raises TypeError: If the provided `gene` is of the wrong type
        """
        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type object was expected, but '{type(gene)}' type object was provided.")
        self.genes.add(gene)
        gene.family = self
        if hasattr(gene, "organism"):
            self._genePerOrg[gene.organism].add(gene)

    def mk_bitarray(self, index: Dict[Organism, int], partition: str = 'all'):
        """Produces a bitarray representing the presence/absence of the family in the pangenome using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.

        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`
        :param partition: partition used to compute bitarray
        """
        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == 'all':
            logging.getLogger("PPanGGOLiN").debug(f"all")
            for org in self.organisms:
                self.bitarray[index[org]] = 1
        elif partition in ['shell', 'cloud']:
            logging.getLogger("PPanGGOLiN").debug(f"shell, cloud")
            if self.named_partition == partition:
                for org in self.organisms:
                    self.bitarray[index[org]] = 1
        elif partition == 'accessory':
            logging.getLogger("PPanGGOLiN").debug(f"accessory")
            if self.named_partition in ['shell', 'cloud']:
                for org in self.organisms:
                    self.bitarray[index[org]] = 1

    def get_org_dict(self) -> Dict[Organism, Set[Gene]]:
        """Returns the organisms and the genes belonging to the gene family

        :return: a dictionnary of organism as key and set of genes as values
        """
        try:
            return self._genePerOrg
        except AttributeError:
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg
        except Exception:
            raise Exception("An unexpected error occurs. Please report in our GitHub")

    def get_genes_per_org(self, org: Organism) -> Set[Gene]:
        """Returns the genes belonging to the gene family in the given Organism

        :param org: Organism to look for

        :return: a set of gene(s)
        """
        try:
            return self._genePerOrg[org]
        except AttributeError:
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg[org]
        except Exception:
            raise Exception("An unexpected error occurs. Please report in our GitHub")

    @property
    def neighbors(self) -> Set[GeneFamily]:
        """Returns all the GeneFamilies that are linked with an edge

        :return: Neighbors
        """
        return set(self._edges.keys())

    @property
    def edges(self) -> List[Edge]:
        """Returns all Edges that are linked to this gene family

        :return: Edges of the gene family
        """
        return list(self._edges.values())

    @property
    def organisms(self) -> Set[Organism]:
        """Returns all the Organisms that have this gene family

        :return: Organisms that have this gene family
        """
        try:
            return set(self._genePerOrg.keys())
        except AttributeError:  # then the genes have been added before they had organisms
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return set(self._genePerOrg.keys())
        except Exception:
            raise Exception("An unexpected error occurs. Please report in our GitHub")
