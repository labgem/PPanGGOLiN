#!/usr/bin/env python3
# coding: utf8

from __future__ import annotations

# installed libraries
import logging
import pdb
from typing import Iterator, Dict
from collections import defaultdict
from bidict import bidict

import gmpy2


class Feature:
    """This is a general class representation of Gene, RNA

    :param identifier: Identifier of the feature given by PPanGGOLiN
    """

    _all_dna = bidict()
    """Feature._all_dna stores all non redondant dna sequences"""

    _featureID2non_redondant_ID = {}
    """Feature._featureID2non_redondant_ID gives the non redondant ID from a feature ID"""

    _non_redondant_ID2featureID_set = defaultdict(set)
    """Feature._non_redondant_ID2featureID_set gives all feature IDs from non a non redondant ID"""

    def __init__(self, identifier: str):
        self.ID = identifier
        self.is_fragment = False
        self.type = ""
        self.start = None
        self.stop = None
        self.strand = None
        self.product = None
        self.name = None
        self.local_identifier = None
        self.organism = None
        self.contig = None
        self._dna = None

    def fill_annotations(self, start: int, stop: int, strand: str, gene_type: str = "", name: str = "",
                         product: str = "", local_identifier: str = ""):
        """
        Fill general annotation for child classes

        :param start: Start position
        :param stop: Stop position
        :param strand: associated strand
        :param gene_type: Type of gene
        :param name: Name of the feature
        :param product: Associated product
        :param local_identifier: Identifier provided by the original file
        """
        self.start = start if isinstance(start, int) else int(start)
        self.stop = stop if isinstance(stop, int) else int(stop)
        self.type = gene_type
        self.strand = strand
        self.product = product
        self.name = name
        self.local_identifier = local_identifier

    def fill_parents(self, organism: Organism, contig: Contig):
        """ Associate object to an organism and a contig

        :param organism: Parent organism
        :param contig: Parent contig
        """
        self.organism = organism
        self.contig = contig

    def add_dna(self, dna):
        """ Add DNA sequence to feature

        :param dna: DNA sequence
        """
        if not isinstance(dna, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(dna)}' type object")
        self._dna = dna

    @property
    def dna(self):
        """ Get dna of a feature

        :return: dna as a str
        """
        if self._dna is not None:
            return (self._dna)
        else:
            return (type(self)._all_dna[type(self)._featureID2non_redondant_ID[self.ID]])

    def add_dna_not_redondant(self, dna):
        """ direclty store dna as not redondant (not thread safe function)

        :param dna: DNA sequence
        """
        self.add_dna(dna)
        self.to_non_redondant_dna()

    def to_non_redondant_dna(self):
        """ compress redondant dna to non redondant dna (not thread safe function)"""
        if self._dna is not None:
            if self._dna in type(self)._all_dna.inverse:
                non_redondant_id = type(self)._all_dna.inverse[self._dna]
            else:
                non_redondant_id = len(type(self)._featureID2non_redondant_ID)
            print(non_redondant_id)
            type(self)._all_dna[non_redondant_id] = self._dna
            type(self)._featureID2non_redondant_ID[self.ID] = non_redondant_id
            type(self)._non_redondant_ID2featureID_set[non_redondant_id].add(self.ID)

class RNA(Feature):
    """Save RNA from genome as an Object with some information for Pangenome

    :param rna_id: Identifier of the rna
    """

    def __init__(self, rna_id: str):
        super().__init__(rna_id)


class Gene(Feature):
    """Save gene from genome as an Object with some information for Pangenome

    :param gene_id: Identifier of the gene
    """
    _all_proteins = bidict()
    """Feature._all_proteins stores all non redondant protein sequences"""

    _proteinIDs2non_redondant_ID = {}
    """Feature._proteinIDs2non_redondant_ID gives the non redondant ID from a protein ID"""

    _non_redondant_ID2proteinIDs_set = defaultdict(set)
    """Feature._non_redondant_ID2proteinIDs_set gives all protein IDs from non a non redondant ID"""

    def __init__(self, gene_id: str):
        super().__init__(gene_id)
        self.position = None
        self.family = None
        self.RGP = set()
        self.genetic_code = None
        self._protein = None

    def __str__(self):
        return str(self.ID)

    def fill_annotations(self, start: int, stop: int, strand: str, gene_type: str = "", name: str = "",
                         product: str = "", local_identifier: str = "", position: int = None, genetic_code: int = 11):
        """
        Fill Gene annotation provide by PPanGGOLiN dependencies

        :param start: Start position
        :param stop: Stop position
        :param strand: associated strand
        :param gene_type: Type of gene
        :param name: Gene name
        :param product: Associated product
        :param local_identifier: Identifier provided by the original file
        :param position: Gene localisation in genome
        :param genetic_code: Genetic code associated to gene
        """
        super().fill_annotations(start, stop, strand, gene_type, name, product, local_identifier)
        self.position = position
        self.genetic_code = genetic_code

    def add_protein(self, protein):
        """ Add animo acid sequence to Gene

        :param protein: protein sequence
        """
        if not isinstance(protein, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(protein)}' type object")
        self._protein = protein

    @property
    def protein(self):
        """ Get protein of a gene

        :return: protein as a str
        """
        if self._protein is not None:
            return (self._protein)
        else:
            return (type(self)._all_proteins[type(self)._proteinIDs2non_redondant_ID[self.ID]])

    def add_protein_not_redondant(self, protein):
        """ direclty store protein as not redondant (not thread safe function)

        :param protein: amino acid sequence
        """
        self.add_protein(protein)
        self.to_non_redondant_protein()

    def to_non_redondant_protein(self):
        """ compress redondant prot to non redondant protein (not thread safe function)"""
        if self._protein is not None:
            if self._protein in type(self)._all_proteins.inverse:
                non_redondant_id = type(self)._all_proteins.inverse[self._protein]
            else:
                non_redondant_id = len(type(self)._proteinIDs2non_redondant_ID)
            print(non_redondant_id)
            type(self)._all_proteins[non_redondant_id] = self._protein
            type(self)._proteinIDs2non_redondant_ID[self.ID] = non_redondant_id
            type(self)._non_redondant_ID2proteinIDs_set[non_redondant_id].add(self.ID)

class Contig:
    """
    Describe the contig content and some information

    :param name: Name of the contig
    :param is_circular: save if the contig is circular
    """
    def __init__(self, name: str, is_circular: bool = False):
        self.name = name
        self.is_circular = is_circular
        self.RNAs = set()  # saving the rna annotations. We're not using them in the vast majority of cases.
        self._genes_start = {}
        self._genes_position = []

    @property
    def genes(self) -> list:
        """ Give the gene content of the contig

        :return: list of gene in contig
        """
        return self._genes_position

    def __str__(self):
        return self.name

    def __iter__(self):
        return iter(self.genes)

    # retrieve gene by start position
    def __getitem__(self, index: int):
        gene = self._genes_start.get(index)
        if not gene:
            if not isinstance(index, int):
                raise TypeError(f"Expected type is int, given type was '{type(index)}'")
            raise IndexError(f"No gene start at the given position {index}")
        return gene

    def add_rna(self, rna: RNA):
        """ Add RNA to contig

        :param rna: RNA object to add
        """
        if not isinstance(rna, RNA):
            raise TypeError(f"'RNA' type was expected but you provided a '{type(rna)}' type object")
        self.RNAs.add(rna)

    def add_gene(self, gene: Gene):
        """ Add gene to Contig

        :param gene: Gene object to add
        """
        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type was expected but you provided a '{type(gene)}' type object")
        if gene.position is None:
            raise TypeError(f"The gene object needs to have its position in the contig filled before adding it")
        while len(self._genes_position) <= gene.position:
            # adding empty values. They should be filled by the end of the parsing.
            # Doing this because genes are not always met in order.
            self._genes_position.append(None)
        self._genes_position[gene.position] = gene
        self._genes_start[gene.start] = gene


class Organism:
    """
    Describe the Genome content and some information

    :param name: Name of the genome
    """
    def __init__(self, name: str):
        self.name = name
        self._contigs_getter = {}
        self.bitarray = None

    @property
    def families(self) -> set:
        """ returns the gene families present in the organism

        :return: set of gene families in organism
        """
        return {gene.family for contig in self.contigs for gene in contig.genes}

    @property
    def genes(self) -> Iterator[Gene]:
        """ Generator to get genes in organism """
        for contig in self.contigs:
            for gene in contig.genes:
                yield gene

    def number_of_genes(self) -> int:
        """ Get number of genes in organism

        :return: Number of gene in organism
        """
        return sum([len(list(contig.genes)) for contig in self.contigs])

    @property
    def contigs(self) -> dict.values:
        """ Get contigs in organism

        :return: values in contig dictionary from organism
        """
        return self._contigs_getter.values()

    def __str__(self):
        return self.name

    def get_contig(self, contig_id: str, is_circular: bool = False):
        """
        Get contig with the given identifier in the organim, if it does not exist in organism,the contig is added

        :param contig_id: Contig idenitifier
        :param is_circular: save if the contig is circular

        :return: the contig with the given identifier
        """
        contig = self._contigs_getter.get(contig_id)
        if contig is None:
            contig = self._create_contig(contig_id, is_circular)
        return contig

    def _create_contig(self, contig_id: str, is_circular: bool = False):
        new_contig = Contig(contig_id, is_circular)
        self._contigs_getter[contig_id] = new_contig
        return new_contig

    def mk_bitarray(self, index: Dict[Organism, int], partition: str = 'all'):
        """Produces a bitarray representing the presence / absence of families in the organism using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.

        :param partition: Filter partition
        :param index: The index computed by :func:`ppanggolin.pan.Pangenome.getIndex`
        """

        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == 'all':
            logging.getLogger().debug(f"all")
            for fam in self.families:
                self.bitarray[index[fam]] = 1
        elif partition in ['shell', 'cloud']:
            logging.getLogger().debug(f"shell, cloud")
            for fam in self.families:
                if fam.named_partition == partition:
                    self.bitarray[index[fam]] = 1
        elif partition == 'accessory':
            logging.getLogger().debug(f"accessory")
            for fam in self.families:
                if fam.named_partition in ['shell', 'cloud']:
                    self.bitarray[index[fam]] = 1
        else:
            raise Exception("There is not any partition corresponding please report a github issue")
