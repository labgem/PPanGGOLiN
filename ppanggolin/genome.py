#!/usr/bin/env python3
# coding: utf8

from __future__ import annotations

# installed libraries
import logging
from typing import Dict, Iterator, Generator

import gmpy2

# local libraries
from ppanggolin.metadata import MetaFeatures


class Feature(MetaFeatures):
    """This is a general class representation of Gene, RNA

    Methods:
    - fill_annotations(): fills general annotation for child classes.
    - fill_parents(): associates the object to an organism and a contig.
    - add_dna(): adds DNA sequence to the feature.

    Fields:
    - ID: Identifier of the feature given by PPanGGOLiN.
    - is_fragment: Boolean value indicating whether the feature is a fragment or not.
    - type: Type of the feature.
    - start: Start position of the feature.
    - stop: Stop position of the feature.
    - strand: Strand associated with the feature.
    - product: Associated product of the feature.
    - name: Name of the feature.
    - local_identifier: Identifier provided by the original file.
    - organism: Parent organism of the feature.
    - contig: Parent contig of the feature.
    - dna: DNA sequence of the feature.
    """
    def __init__(self, identifier: str):
        """Constructor Method

        :param identifier: identifier of the feature
        """
        assert isinstance(identifier, str), "Expected identifier should be a string"
        if identifier == '':
            raise ValueError("Identifier should not be empty")
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
        self._organism = None
        self._contig = None
        self.dna = None

    @property
    def organism(self) -> Organism:
        """Return organism that Feature belongs to.

        :return: Organism of the feature
        """
        return self._organism

    @organism.setter
    def organism(self, organism: Organism):
        if not isinstance(organism, Organism):
            raise TypeError(f'Expected type Organism, got {type(organism)}')
        self._organism = organism

    @property
    def contig(self) -> Contig:
        """Return contig that Feature belongs to.

        :return: Contig of the feature
        """
        return self._contig

    @contig.setter
    def contig(self, contig: Contig):
        if not isinstance(contig, Contig):
            raise TypeError(f'Expected type Contig, got {type(contig)}')
        self._contig = contig

    @property
    def length(self) -> int:
        """Return gene length

        :return: gene length
        """
        if self.start is not None:
            if self.stop is not None:
                return self.stop - self.start
            else:
                raise ValueError("Stop is not known")
        else:
            raise ValueError("Start is not known")

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

        :raises TypeError: If attribute value not correspond to expected type
        :raises ValueError: If strand is not '+' or '-'
        """
        if not isinstance(start, int):
            raise TypeError("Start should be int")
        if not isinstance(stop, int):
            raise TypeError("Stop should be int")
        if not isinstance(strand, str):
            raise TypeError("Strand should be str")
        if not isinstance(gene_type, str):
            raise TypeError("Gene type should be str")
        if not isinstance(name, str):
            raise TypeError("Name should be str")
        if not isinstance(product, str):
            raise TypeError("Product should be str")
        if not isinstance(local_identifier, str):
            raise TypeError("Local identifier should be str")
        if strand not in ["+", "-"]:
            raise ValueError("Strand should be + or -")
        self.start = start
        self.stop = stop
        self.strand = strand
        self.type = gene_type
        self.product = product
        self.name = name
        self.local_identifier = local_identifier

    def fill_parents(self, organism: Organism = None, contig: Contig = None):
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
        assert isinstance(dna, str), f"'str' type was expected but you provided a '{type(dna)}' type object"
        self.dna = dna


class RNA(Feature):
    """Save RNA from genome as an Object with some information for Pangenome

    :param rna_id: Identifier of the rna
    """

    def __init__(self, rna_id: str):
        super().__init__(rna_id)


class Gene(Feature):
    """Save gene from genome as an Object with some information for Pangenome

    Methods:
    - fill_annotations(): fills general annotation for the gene object and adds additional attributes such as
    position and genetic code.
    - add_protein(): adds the protein sequence corresponding to the translated gene to the object.

    Fields:
    - position: the position of the gene in the genome.
    - family: the family that the gene belongs to.
    - RGP: a set of resistance gene profiles associated with the gene.
    - genetic_code: the genetic code associated with the gene.
    - protein: the protein sequence corresponding to the translated gene.
    """
    def __init__(self, gene_id: str):
        """Constructor method

        :param gene_id: Identifier of the gene
        """
        super().__init__(gene_id)
        self.position = None
        self._family = None
        self._RGP = None
        self.genetic_code = None
        self.protein = None

    def __str__(self) -> str:
        return str(self.ID)

    @property
    def family(self):
        """Return GeneFamily that Gene belongs to.

        :return: Gene family of the gene
        :rtype: GeneFamily
        """
        return self._family

    @family.setter
    def family(self, family):
        from ppanggolin.geneFamily import GeneFamily
        if not isinstance(family, GeneFamily):
            raise TypeError(f'Expected type Organism, got {type(family)}')
        self._family = family

    @property
    def RGP(self):
        """Return the RGP that gene belongs to

        :return: RGP fo the Gene
        :rtype: Region
        """
        return self._RGP

    @RGP.setter
    def RGP(self, RGP):
        from ppanggolin.region import Region
        if not isinstance(RGP, Region):
            raise TypeError(f'Expected type Organism, got {type(RGP)}')
        self._RGP = RGP

    def fill_annotations(self, position: int = None, genetic_code: int = 11, **kwargs):
        """Fill Gene annotation provide by PPanGGOLiN dependencies

        :param position: Gene localisation in genome
        :param genetic_code: Genetic code associated to gene
        :param kwargs: look at Feature.fill_annotations methods
        """
        super().fill_annotations(**kwargs)
        if position is not None and not isinstance(position, int):
            raise TypeError("position should be an integer")
        if not isinstance(genetic_code, int):
            raise TypeError("Genetic code should be an integer")
        self.position = position
        self.genetic_code = genetic_code

    def add_protein(self, protein: str):
        """Add protein sequence corresponding to translated gene

        :param protein: Protein sequence

        :raise TypeError: Protein sequence must be a string
        """
        if not isinstance(protein, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(protein)}' type object")
        self.protein = protein


class Contig:
    """
    Describe the contig content and some information
    Methods:
    - genes(self) -> list: Returns a list of gene objects present in the contig.
    - add_rna(self, rna: RNA): Adds an RNA object to the contig.
    - add_gene(self, gene: Gene): Adds a gene object to the contig.

    Fields:
    - name: Name of the contig.
    - is_circular: Boolean value indicating whether the contig is circular or not.
    - RNAs: Set of RNA annotations present in the contig.
    """

    def __init__(self, name: str, is_circular: bool = False):
        """Constructor method
        :param name: Name of the contig
        :param is_circular: save if the contig is circular
        """
        self.name = name
        self.is_circular = is_circular
        self._rnaGetter = set()  # saving the rna annotations. We're not using them in the vast majority of cases.
        self._genesGetter = {}
        self._genes_position = []

    def __str__(self) -> str:
        return self.name

    def __len__(self):
        return len(self._genes_position)

    def __setitem__(self, start: int, gene: Gene):
        """ Set gene to Contig

        :param start: Start position of the gene
        :param gene: Gene object to add
        """
        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type was expected but you provided a '{type(gene)}' type object")
        if start in self._genesGetter:
            raise ValueError(f"Gene with start position {start} already exists in the contig")
        if gene.position is None:
            raise AttributeError("The gene object needs to have its position in the contig filled before adding it")
        # adding empty values. They should be filled by the end of the parsing.
        # Doing this because genes are not always met in order.
        self._genes_position.extend([None] * (gene.position - len(self._genes_position) + 1))
        self._genes_position[gene.position] = gene
        self._genesGetter[gene.start] = gene

    # retrieve gene by start position
    def __getitem__(self, index: int) -> Gene:
        if not isinstance(index, int):
            raise TypeError(f"Expected type is int, given type was '{type(index)}'")
        return self._genes_position[index]

    def get_genes(self, begin: int, end: int):
        if not isinstance(begin, int) or not isinstance(end, int):
            raise TypeError(f"Expected type is int, given type was '{type(begin)}, {type(end)}'")
        if end < begin:
            raise ValueError("End position is lower than begin position")
        else:
            return self._genes_position[begin: end]

    @property
    def genes(self) -> list:
        """ Give the gene content of the contig

        :return: list of gene in contig
        """
        for gene in self._genes_position:
            yield gene

    @property
    def RNAs(self) -> Generator[RNA, None, None]:
        """Return all the RNA in the contig

        :return: Generator of RNA
        """
        for rna in self._rnaGetter:
            yield rna

    def add_rna(self, rna: RNA):
        """ Add RNA to contig

        :param rna: RNA object to add
        """
        if not isinstance(rna, RNA):
            raise TypeError(f"'RNA' type was expected but you provided a '{type(rna)}' type object")
        if rna in self._rnaGetter:
            raise KeyError(f"RNA with the id: {rna.ID} already exist in contig {self.name}")
        self._rnaGetter.add(rna)


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
        return sum([len(contig) for contig in self.contigs])

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
