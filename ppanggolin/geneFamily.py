#!/usr/bin/env python3
# coding: utf8

# default libraries
from __future__ import annotations
from collections import defaultdict
import logging

# installed libraries
from typing import Dict, Generator, List, Set, Tuple, Union

import gmpy2

# local libraries
from ppanggolin.edge import Edge
from ppanggolin.genome import Gene, Organism
from ppanggolin.metadata import Metadata


class GeneFamily:
    """
    This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.

    :param family_id: The internal identifier to give to the gene family
    :type family_id: any
    :param name: The name of the gene family (to be printed in output files)
    :type name: str
    """

    def __init__(self, family_id: int, name: str):
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
        self._metadataGetter = {}  # Key = source, Value = ordered list of the best annotation for one source

    def add_sequence(self, seq: str):
        """Assigns a protein sequence to the gene family.

        :param seq: the sequence to add to the gene family
        """
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

        :param index: The index computed by :func:`ppanggolin.pan.Pangenome.getIndex`
        :param partition: partition used to compute bitarray
        """
        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == 'all':
            logging.getLogger().debug("all")
            for org in self.organisms:
                self.bitarray[index[org]] = 1
        elif partition in ['shell', 'cloud']:
            logging.getLogger().debug("shell, cloud")
            if self.named_partition == partition:
                for org in self.organisms:
                    self.bitarray[index[org]] = 1
        elif partition == 'accessory':
            logging.getLogger().debug("accessory")
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

    @property
    def metadata(self) -> Generator[Metadata, None, None]:
        """Generate annotations in gene families

        :return: Gene family annotation"""
        for annot_list in self._metadataGetter.values():
            for annotation in annot_list:
                yield annotation

    @property
    def sources(self) -> List[str]:
        """ Get all metadata source in gene family

        :return: List of annotation source
        """
        return list(self._metadataGetter.keys())

    def max_metadata_by_source(self) -> Tuple[str, int]:
        """Get the maximum number of annotation for one source

        :return: Name of the source with the maximum annotation and the number of annotation corresponding
        """
        max_annot = 0
        max_source = None
        for source, annotations in self._metadataGetter.items():
            if len(annotations) > max_annot:
                max_annot = len(annotations)
                max_source = source
        return max_source, max_annot

    def get_source(self, name: str) -> Union[List[Metadata], None]:
        """ Get the annotation for a specific source in gene family

        :param name: Name of the source

        :return: All the annotation from the source if exist else None
        """
        return self._metadataGetter[name] if name in self.sources else None

    def get_metadata(self, value: Union[str, int, float, List[str], List[int], List[float]],
                     accession: Union[List[str], str]) -> Generator[Metadata, None, None]:
        """Get annotation by name or accession in gene family

        :param value: Names of annotation searched
        :param accession: Accession number of annotation searched

        :return: annotation searched
        """
        assert value is not None and accession is not None
        value = value if isinstance(value, list) else [value]
        accession = accession if isinstance(accession, list) else [accession]

        for annotation in self.metadata:
            if annotation.value in value or annotation.accession in accession:
                yield annotation

    def add_metadata(self, source: str, metadata: Metadata):
        """ Add annotation to gene family

        :param source: Name of database source
        :param metadata: Identifier of the annotation
        """
        source_annot = self.get_source(source)
        same_value = False
        if source_annot is not None:
            index_annot = 0
            insert_bool = False
            while index_annot < len(source_annot):
                current_annot = source_annot[index_annot]
                if current_annot.value == metadata.value:
                    same_value = True
                if current_annot.score is not None and metadata.score is not None:
                    if current_annot.score < metadata.score:
                        if same_value:
                            source_annot[index_annot] = metadata
                        else:
                            source_annot.insert(index_annot, metadata)
                            insert_bool = True
                    elif current_annot.score == metadata.score:
                        if current_annot.e_val is not None and metadata.e_val is not None:
                            if current_annot.e_val > metadata.e_val:
                                if same_value:
                                    source_annot[index_annot] = metadata
                                else:
                                    source_annot.insert(index_annot, metadata)
                                    insert_bool = True
                elif current_annot.e_val is not None and metadata.e_val is not None:
                    if current_annot.e_val > metadata.e_val:
                        if same_value:
                            source_annot[index_annot] = metadata
                        else:
                            source_annot.insert(index_annot, metadata)
                            insert_bool = True
                if not insert_bool and not same_value:
                    index_annot += 1
                else:
                    break
            if not insert_bool and not same_value:
                source_annot.append(metadata)
        else:
            self._metadataGetter[source] = [metadata]
