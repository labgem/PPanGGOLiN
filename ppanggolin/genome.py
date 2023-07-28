#!/usr/bin/env python3
# coding: utf8

from __future__ import annotations

# installed libraries
import logging
from typing import Dict, Iterator

import gmpy2

# local libraries
from ppanggolin.metadata import MetaFeatures


class Feature(MetaFeatures):
    """This is a general class representation of Gene, RNA

    :param identifier: Identifier of the feature given by PPanGGOLiN
    """

    def __init__(self, identifier: str):
        super().__init__()
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
        self.dna = None

    @property
    def length(self) -> int:
        """Return gene length

        :return: gene length
        """
        return self.stop - self.start

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

        :raise TypeError: DNA sequence must be a string
        """
        if not isinstance(dna, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(dna)}' type object")
        self.dna = dna


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

    def __init__(self, gene_id: str):
        super().__init__(gene_id)
        self.position = None
        self.family = None
        self.RGP = set()
        self.genetic_code = None
        self.protein = None

    def __str__(self) -> str:
        return str(self.ID)

    def fill_annotations(self, start: int, stop: int, strand: str, gene_type: str = "", name: str = "",
                         product: str = "", local_identifier: str = "", position: int = None, genetic_code: int = 11):
        """
        Fill Gene annotation provide by PPanGGOLiN dependencies

        :param start: Start position
        :param stop: Stop position
        :param strand: associated strand
        :param gene_type: Type of the gene
        :param name: Gene name
        :param product: Associated product
        :param local_identifier: Identifier provided by the original file
        :param position: Gene localisation in genome
        :param genetic_code: Genetic code associated to gene
        """
        super().fill_annotations(start, stop, strand, gene_type, name, product, local_identifier)
        self.position = position
        self.genetic_code = genetic_code

    def add_protein(self, protein: str):
        """ Add protein sequence corresponding to translated gene

        :param protein: Protein sequence

        :raise TypeError: Protein sequence must be a string
        """
        if not isinstance(protein, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(protein)}' type object")
        self.protein = protein


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

    def __str__(self) -> str:
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
            raise TypeError("The gene object needs to have its position in the contig filled before adding it")
        while len(self._genes_position) <= gene.position:
            # adding empty values. They should be filled by the end of the parsing.
            # Doing this because genes are not always met in order.
            self._genes_position.append(None)
        self._genes_position[gene.position] = gene
        self._genes_start[gene.start] = gene


class Organism(MetaFeatures):
    """
    Describe the Genome content and some information

    :param name: Name of the genome
    """
    def __init__(self, name: str):
        super().__init__()
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
        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`
        """

        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == 'all':
            logging.getLogger("PPanGGOLiN").debug(f"all")
            for fam in self.families:
                self.bitarray[index[fam]] = 1
        elif partition in ['shell', 'cloud']:
            logging.getLogger("PPanGGOLiN").debug(f"shell, cloud")
            for fam in self.families:
                if fam.named_partition == partition:
                    self.bitarray[index[fam]] = 1
        elif partition == 'accessory':
            logging.getLogger("PPanGGOLiN").debug(f"accessory")
            for fam in self.families:
                if fam.named_partition in ['shell', 'cloud']:
                    self.bitarray[index[fam]] = 1
        else:
            raise Exception("There is not any partition corresponding please report a github issue")
