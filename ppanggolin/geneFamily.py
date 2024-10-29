#!/usr/bin/env python3

# default libraries
from __future__ import annotations
from collections import defaultdict
import logging

# installed libraries
from typing import Dict, Generator, Set
import gmpy2

# local libraries
from ppanggolin.edge import Edge
from ppanggolin.genome import Gene, Organism
from ppanggolin.metadata import MetaFeatures


class GeneFamily(MetaFeatures):
    """
    This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.

    Methods:
        - named_partition: returns a meaningful name for the partition associated with the family.
        - neighbors: returns all the GeneFamilies that are linked with an edge.
        - edges: returns all Edges that are linked to this gene family.
        - genes: returns all the genes associated with the family.
        - organisms: returns all the Organisms that have this gene family.
        - spots: returns all the spots associated with the family.
        - modules: returns all the modules associated with the family.
        - number_of_neighbor: returns the number of neighbor GeneFamilies.
        - number_of_edges: returns the number of edges.
        - number_of_genes: returns the number of genes.
        - number_of_organisms: returns the number of organisms.
        - number_of_spots: returns the number of spots.
        - set_edge: sets an edge between the current family and a target family.
        - add_sequence: assigns a protein sequence to the gene family.
        - add_gene: adds a gene to the gene family and sets the gene's family accordingly.
        - add_spot: adds a spot to the gene family.
        - add_module: adds a module to the gene family.
        - Mk_bitarray: produces a bitarray representing the presence/absence of the family in the pangenome using the provided index.
        - get_org_dict: returns a dictionary of organisms as keys and sets of genes as values.
        - get_genes_per_org: returns the genes belonging to the gene family in the given organism.

    Fields:
        - name: the name of the gene family.
        - ID: the internal identifier of the gene family.
        - removed: a boolean indicating whether the family has been removed from the main graph.
        - sequence: the protein sequence associated with the family.
        - Partition: the partition associated with the family.
    """

    def __init__(self, family_id: int, name: str):
        # TODO edges as genes in contig to get and set
        """Constructor method
        :param family_id: The internal identifier to give to the gene family
        :type family_id: any
        :param name: The name of the gene family (to be printed in output files)
        :type name: str
        """
        assert isinstance(family_id, int), "GeneFamily object id should be an integer"
        assert isinstance(name, str), "GeneFamily object name should be a string"
        assert name != "", "GeneFamily object cannot be created with an empty name"

        super().__init__()
        self.name = str(name)
        self.ID = family_id
        self._representative = None
        self._edges_getter = {}
        self._genePerOrg = defaultdict(set)
        self._genes_getter = {}
        self.removed = False  # for the repeated family not added in the main graph
        self.sequence = ""
        self._partition = None
        self._spots = set()
        self._module = None
        self.bitarray = None

    def __repr__(self) -> str:
        """Family representation"""
        return f"{self.ID}: {self.name}"

    def __len__(self) -> int:
        """Get the number of genes in the family

        :return: The length of the list of genes
        """
        return len(self._genes_getter)

    def __setitem__(self, identifier: str, gene: Gene):
        """Set gene to Gene Family

        :param identifier: ID of the gene
        :param gene: Gene object to add

        :raises TypeError: If the gene is not instance Gene
        :raises TypeError: If the identifier is not instance string
        :raises ValueError: If a gene in getter already exists at the name
        """
        # TODO look at change start for position

        if not isinstance(gene, Gene):
            raise TypeError(
                f"'Gene' type was expected but you provided a '{type(gene)}' type object"
            )
        if not isinstance(identifier, str):
            raise TypeError(
                f"Gene ID should be a string. You provided a '{type(identifier)}' type object"
            )
        if identifier in self._genes_getter:
            raise KeyError(
                f"Gene with name {identifier} already exists in the gene family"
            )
        self._genes_getter[identifier] = gene

    # TODO define eq function

    # retrieve gene by start position
    def __getitem__(self, identifier: str) -> Gene:
        """Get the gene for the given name

        :param identifier: ID of the gene in the gene family

        :return:  Wanted gene

        :raises TypeError: If the identifier is not instance string
        :raises KeyError: Gene with the given identifier does not exist in the contig
        """
        if not isinstance(identifier, str):
            raise TypeError(
                f"Gene ID should be a string. You provided a '{type(identifier)}' type object"
            )
        try:
            gene = self._genes_getter[identifier]
        except KeyError:
            raise KeyError(
                f"Gene with the ID: {identifier} does not exist in the family"
            )
        else:
            return gene

    def __delitem__(self, identifier: str):
        """Remove the gene for the given name in the gene family

        :param position: ID of the gene in the family

        :raises TypeError: If the identifier is not instance string
        :raises KeyError: Gene with the given identifier does not exist in the contig
        """
        if not isinstance(identifier, str):
            raise TypeError(
                f"Gene ID should be a string. You provided a '{type(identifier)}' type object"
            )
        try:
            del self._genes_getter[identifier]
        except KeyError:
            raise KeyError(
                f"Gene with the name: {identifier} does not exist in the family"
            )

    def add(self, gene: Gene):
        """Add a gene to the gene family, and sets the gene's :attr:family accordingly.

        :param gene: The gene to add

        :raises TypeError: If the provided `gene` is of the wrong type
        """
        if not isinstance(gene, Gene):
            raise TypeError(
                f"'Gene' type object was expected, but '{type(gene)}' type object was provided."
            )
        self[gene.ID] = gene
        gene.family = self
        if gene.organism is not None and gene.organism in self._genePerOrg:
            # TODO try to remove the second condition and check if projection is working
            self._genePerOrg[gene.organism].add(gene)

    def get(self, identifier: str) -> Gene:
        """Get a gene by its name

        :param identifier: ID of the gene

        :return: Wanted gene

        :raises TypeError: If the identifier is not instance string
        """
        if not isinstance(identifier, str):
            raise TypeError(
                f"Gene ID should be a string. You provided a '{type(identifier)}' type object"
            )
        return self[identifier]

    def remove(self, identifier):
        """Remove a gene by its name

        :param identifier: Name of the gene

        :return: Wanted gene

        :raises TypeError: If the identifier is not instance string
        """
        if not isinstance(identifier, str):
            raise TypeError(
                f"Gene ID should be a string. You provided a '{type(identifier)}' type object"
            )
        del self[identifier]

    @property
    def representative(self) -> Gene:
        """Get the representative gene of the family

        :return: The representative gene of the family
        """
        if self._representative is None:
            raise Exception("Representative gene has not been set")
        return self._representative

    @representative.setter
    def representative(self, gene: Gene) -> None:
        """Set the representative gene of the family"""
        if not isinstance(gene, Gene):
            raise TypeError(
                f"Representative gene should be a Gene. Found a '{type(gene)}' type object"
            )
        self._representative = gene

    def contains_gene_id(self, identifier):
        """
        Check if the family contains already a gene id

        :param identifier: ID of the gene

        :return: True if it contains False if it does not

        :raises TypeError: If the identifier is not instance string
        """
        if not isinstance(identifier, str):
            raise TypeError(
                f"Gene ID should be a string. You provided a '{type(identifier)}' type object"
            )

        return identifier in self._genes_getter

    # TODO define __eq__
    @property
    def partition(self):
        return self._partition if self._partition is not None else ""

    @partition.setter
    def partition(self, partition: str):
        self._partition = partition

    @property
    def named_partition(self) -> str:
        """Reads the partition attribute and returns a meaningful name

        :return: The partition name of the gene family

        :raises ValueError: If the gene family has no partition assigned
        """
        if self.partition == "":
            raise ValueError("The gene family has not been associated to a partition.")
        if self.partition.startswith("P"):
            return "persistent"
        elif self.partition.startswith("C"):
            return "cloud"
        elif self.partition.startswith("S"):
            return "shell"
        else:
            return "undefined"

    @property
    def edges(self) -> Generator[Edge, None, None]:
        """Returns all Edges that are linked to this gene family

        :return: Edges of the gene family
        """
        yield from self._edges_getter.values()

    @property
    def neighbors(self) -> Generator[GeneFamily, None, None]:
        """Returns all the GeneFamilies that are linked with an edge

        :return: Neighbors
        """
        yield from self._edges_getter.keys()

    @property
    def genes(self):
        """Return all the genes belonging to the family

        :return: Generator of genes
        """
        yield from self._genes_getter.values()

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """Returns all the Organisms that have this gene family

        :return: Organisms that have this gene family
        """
        if len(self._genePerOrg) == 0:
            _ = self.get_org_dict()
        yield from self._genePerOrg.keys()

    @property
    def spots(self) -> Generator[Spot, None, None]:
        """Return all the spots belonging to the family

        :return: Generator of spots
        """
        yield from self._spots

    @property
    def module(self) -> Module:
        """Return all the modules belonging to the family

        :return: Generator of modules
        """
        return self._module

    @module.setter
    def module(self, module: Module):
        """Set the modules belonging to the family
        :param module: module to set
        """
        if not isinstance(module, Module):
            raise TypeError("Module object is expected to object of Module class")
        self._module = module

    @property
    def number_of_neighbors(self) -> int:
        """Get the number of neighbor for the current gene family"""
        return len(self._edges_getter.keys())

    @property
    def number_of_edges(self) -> int:
        """Get the number of edges for the current gene family"""
        return len(self._edges_getter.values())

    @property
    def number_of_genes(self) -> int:
        """Get the number of genes for the current gene family"""
        return len(self._genes_getter)

    @property
    def number_of_organisms(self) -> int:
        """Get the number of organisms for the current gene family"""
        if len(self._genePerOrg) == 0:
            _ = self.get_org_dict()
        return len(self._genePerOrg.keys())

    @property
    def number_of_spots(self) -> int:
        """Get the number of spots for the current gene family"""
        return len(self._spots)

    @property
    def has_module(self) -> bool:
        """
        Check if the family is in a module

        return True if it has a module else False
        """
        return self._module is not None

    def set_edge(self, target: GeneFamily, edge: Edge):
        """Set the edge between the gene family and another one

        :param target: Neighbor family
        :param edge: Edge connecting families
        """
        self._edges_getter[target] = edge

    def get_edge(self, target: GeneFamily) -> Edge:
        """Get the edge by the target gene family neighbor"""
        return self._edges_getter[target]

    def add_sequence(self, seq: str):
        """Assigns a protein sequence to the gene family.

        :param seq: The sequence to add to the gene family
        """
        assert isinstance(seq, str), "Sequence must be a string"

        self.sequence = seq

    def add_spot(self, spot: Spot):
        """Add the given spot to the family

        :param spot: Spot belonging to the family
        """
        from ppanggolin.region import Spot  # prevent circular import error

        if not isinstance(spot, Spot):
            raise TypeError(f"A spot object is expected, you give a {type(spot)}")
        self._spots.add(spot)

    def set_module(self, module: Module):
        """Add the given module to the family

        :param module: Module belonging to the family
        """
        from ppanggolin.region import Module  # prevent circular import error

        if not isinstance(module, Module):
            raise TypeError(f"A module object is expected, you give a {type(module)}")
        self._module = module

    def mk_bitarray(self, index: Dict[Organism, int], partition: str = "all"):
        """Produces a bitarray representing the presence/absence of the family in the pangenome using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.

        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`
        :param partition: partition used to compute bitarray
        """
        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == "all":
            logging.getLogger("PPanGGOLiN").debug("all")
            for org in self.organisms:
                self.bitarray[index[org]] = 1
        elif partition in ["shell", "cloud"]:
            logging.getLogger("PPanGGOLiN").debug("shell, cloud")
            if self.named_partition == partition:
                for org in self.organisms:
                    self.bitarray[index[org]] = 1
        elif partition == "accessory":
            logging.getLogger("PPanGGOLiN").debug("accessory")
            if self.named_partition in ["shell", "cloud"]:
                for org in self.organisms:
                    self.bitarray[index[org]] = 1

    def get_org_dict(self) -> Dict[Organism, Set[Gene]]:
        """Returns the organisms and the genes belonging to the gene family

        :return: A dictionary of organism as key and set of genes as values
        """
        if len(self._genePerOrg) == 0:
            for gene in self.genes:
                if gene.organism is None:
                    raise AttributeError(f"Gene: {gene.name} is not fill with genome")
                self._genePerOrg[gene.organism].add(gene)
        return self._genePerOrg

    def get_genes_per_org(self, org: Organism) -> Generator[Gene, None, None]:
        """Returns the genes belonging to the gene family in the given Organism

        :param org: Organism to look for

        :return: A set of gene(s)
        """
        if len(self._genePerOrg) == 0:
            _ = self.get_org_dict()
        if org not in self._genePerOrg:
            raise KeyError(
                f"Genome {org.name} does not have the gene family: {self.name}"
            )
        yield from self._genePerOrg[org]

    def is_single_copy(self, dup_margin: float, exclude_fragment: bool) -> bool:
        """
        Checks if the gene family is considered single copy based on the provided criteria.

        :param dup_margin: The maximum allowed duplication margin for a gene family to be considered single copy.
        :param exclude_fragment: A boolean indicating whether to exclude fragments when determining single copy families.
        :return: A boolean indicating whether the gene family is single copy.
        """

        return self.duplication_ratio(exclude_fragment) < dup_margin

    def duplication_ratio(self, exclude_fragment: bool) -> bool:
        """
        Checks if the gene family is considered single copy based on the provided criteria.

        :param dup_margin: The maximum allowed duplication margin for a gene family to be considered single copy.
        :param exclude_fragment: A boolean indicating whether to exclude fragments when determining single copy families.
        :return: A boolean indicating whether the gene family is single copy.
        """
        orgs_with_fam_in_multicopy = 0

        # Check if the family is in multicopy in all organisms
        for fam_genes_in_org in self.get_org_dict().values():
            if exclude_fragment:
                genes_count = len(
                    [gene for gene in fam_genes_in_org if not gene.is_fragment]
                )
            else:
                genes_count = len(fam_genes_in_org)

            if genes_count > 1:
                orgs_with_fam_in_multicopy += 1

        return orgs_with_fam_in_multicopy / self.number_of_organisms
