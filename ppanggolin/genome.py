#!/usr/bin/env python3
# coding: utf8

from __future__ import annotations

# installed libraries
import logging
from typing import Dict, Generator, List

import gmpy2

# local libraries
from ppanggolin.metadata import MetaFeatures


class Feature(MetaFeatures):
    """This is a general class representation of Gene, RNA

    Methods:
    - fill_annotations: fills general annotation for child classes.
    - fill_parents: associates the object to an organism and a contig.
    - Add_sequence: adds a sequence to the feature.

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

        :param identifier: Identifier of the feature
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

    def __str__(self) -> str:
        return str(self.ID)

    def __len__(self) -> int:
        """Return gene length

        :return: gene length

        :raises ValueError: If start or stop are not defined in gene
        """
        if self.start is not None:
            if self.stop is not None:
                return self.stop - self.start + 1
            else:
                raise ValueError("Stop is not known")
        else:
            raise ValueError("Start is not known")

    @property
    def organism(self) -> Organism:
        """Return organism that Feature belongs to.

        :return: Organism of the feature
        """
        return self._organism

    @organism.setter
    def organism(self, organism: Organism):
        """Set the organism to the Feature

        :param organism: Organism belonging to the feature
        """
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
        """Set the contig to the Feature

        :param contig: Contig linked to the feature
        """
        if not isinstance(contig, Contig):
            raise TypeError(f'Expected type Contig, got {type(contig)}')
        self._contig = contig

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

        :raises TypeError: If attribute value does not correspond to the expected type
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
        if organism is not None:
            # TODO test type
            self.organism = organism
            if contig is not None:
                self.contig = contig
        else:
            if contig is not None:
                self.contig = contig
            else:
                raise AssertionError("You should provide at least organism or contig")

    def add_sequence(self, sequence):
        """Add a sequence to feature

        :param sequence: Sequence corresponding to the feature

        :raise AssertionError: Sequence must be a string
        """
        assert isinstance(sequence, str), f"'str' type was expected but you provided a '{type(sequence)}' type object"
        self.dna = sequence


class RNA(Feature):
    """Save RNA from genome as an Object with some information for Pangenome

    :param rna_id: Identifier of the rna
    """

    def __init__(self, rna_id: str):
        super().__init__(rna_id)


class Gene(Feature):
    """Save gene from the genome as an Object with some information for Pangenome

    Methods:
    - fill_annotations: fills general annotation for the gene object and adds additional attributes such as
    position and genetic code.
    - Add_protein: adds the protein sequence corresponding to the translated gene to the object.

    Fields:
    - position: the position of the gene in the genome.
    - family: the family that the gene belongs to.
    - RGP: A putative Region of Plasticity that contains the gene. 
    - genetic_code: the genetic code associated with the gene.
    - Protein: the protein sequence corresponding to the translated gene.
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

    @property
    def family(self):
        """Return GeneFamily that Gene belongs to.

        :return: Gene family of the gene
        :rtype: GeneFamily
        """
        return self._family

    @family.setter
    def family(self, family):
        """Set the GeneFamily blonging to the gene

        :param family: Gene family linked to the gene
        """
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
    def RGP(self, region):
        """Set the Region blonging to the gene

        :param region: Region linked to the gene
        """
        from ppanggolin.region import Region
        if not isinstance(region, Region):
            raise TypeError(f'Expected type Organism, got {type(region)}')
        self._RGP = region

    def fill_annotations(self, position: int = None, genetic_code: int = 11, **kwargs):
        """Fill Gene annotation provide by PPanGGOLiN dependencies

        :param position: Gene localization in genome
        :param genetic_code: Genetic code associated to gene
        :param kwargs: look at Feature.fill_annotations methods

        :raises TypeError: If position or genetic code value is not instance integers
        """
        super().fill_annotations(**kwargs)
        if position is not None and not isinstance(position, int):
            raise TypeError("position should be an integer")
        if not isinstance(genetic_code, int):
            raise TypeError("Genetic code should be an integer")
        self.position = position
        self.genetic_code = genetic_code

    def add_protein(self, protein: str):
        """Add a protein sequence corresponding to translated gene

        :param protein: Protein sequence

        :raise TypeError: Protein sequence must be a string
        """
        if not isinstance(protein, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(protein)}' type object")
        self.protein = protein


class Contig(MetaFeatures):
    """
    Describe the contig content and some information
    Methods:
    - genes: Returns a list of gene objects present in the contig.
    - add_rna: Adds an RNA object to the contig.
    - add_gene: Adds a gene object to the contig.

    Fields:
    - name: Name of the contig.
    - is_circular: Boolean value indicating whether the contig is circular or not.
    - RNAs: Set of RNA annotations present in the contig.
    """

    def __init__(self, identifier: int, name: str, is_circular: bool = False):
        """Constructor method

        :param name: Name of the contig
        :param is_circular: saves if the contig is circular
        """
        super().__init__()
        self.ID = identifier
        self.name = name
        self.is_circular = is_circular
        self._rna_getter = set()  # Saving the rna annotations. We're not using them in the vast majority of cases.
        self._genes_getter = {}
        self._genes_position = []
        self._organism = None
        self._length = None

    def __str__(self) -> str:
        return self.name

    def __setitem__(self, start: int, gene: Gene):
        """ Set gene to Contig

        :param start: Start position of the gene
        :param gene: Gene object to add

        :raises TypeError: If the gene is not instance Gene
        :raises ValueError: If a gene in getter already exists at the start
        :raises AttributeError: If the gene position in the contig is not fill
        """
        # TODO look at change start for position

        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type was expected but you provided a '{type(gene)}' type object")
        if start in self._genes_getter:
            raise ValueError(f"Gene '{self._genes_getter[start].ID}' with start position {start} already exists in the "
                             f"contig '{self.name}' {f'from organism {self.organism}' if self.organism is not None else ''}, "
                             f"cannot add gene '{gene.ID}' {f'from organism {gene.organism}' if gene.organism is not None else ''}")
        if gene.position is None:
            raise AttributeError("The gene object needs to have its position in the contig filled before adding it")
        # Adding empty values.
        # They should be filled by the end of the parsing.
        # Doing this because genes are not always met in order.
        self._genes_position.extend([None] * (gene.position - len(self._genes_position) + 1))
        self._genes_position[gene.position] = gene
        self._genes_getter[gene.start] = gene

    # TODO define eq function

    @property
    def length(self):
        if self._length is None:
            logging.getLogger("PPanGGOLiN").warning("Contig length is unknown")
        return self._length

    @length.setter
    def length(self, contig_len: int):
        if not isinstance(contig_len, int):
            raise TypeError("Contig length is expected to be an integer")
        if contig_len < 0:
            raise ValueError("Contig length must be positive")

        if self._length is None:
            self._length = contig_len
        elif self.length != contig_len:
            logging.getLogger("PPanGGOLiN").debug(f"Known contig length = {self.length}, new length = {contig_len}")
            raise ValueError('Attempting to define a contig length different from the previously defined value.')
        

    def __len__(self):
        return self.length

    # retrieve gene by start position
    def __getitem__(self, position: int) -> Gene:
        """Get the gene for the given position

        :param position: Position of the gene in the contig

        :return:  Wanted gene for the position

        :raises TypeError: If position is not an integer
        """
        if not isinstance(position, int):
            raise TypeError(f"Expected type is int, given type was '{type(position)}'")
        try:
            return self._genes_position[position]
        except KeyError:
            raise KeyError("Position of the gene in the contig does not exist")

    def __delitem__(self, position):
        """Remove the gene for the given position in the contig

        :param position: Position of the gene in the contig

        :raises KeyError: Gene at the given position does not exist in the contig
        """
        if not isinstance(position, int):
            raise TypeError(f"Expected type is int, given type was '{type(position)}'")
        try:
            del self._genes_position[position]
        except KeyError:
            raise KeyError("Position of the gene in the contig does not exist")

    def add(self, gene: Gene):
        """Add a gene to the contig

        :param gene: Gene to add

        :raises TypeError: Region is not an instance Region
        """
        if not isinstance(gene, Gene):
            raise TypeError(f"Unexpected class / type for {type(gene)} when adding it to a contig")
        if gene.start is None:
            raise AttributeError(f'Gene {gene.name} is not fill with start')
        if gene.position is None:
            raise AttributeError(f'Gene {gene.name} is not fill with position')
        self[gene.start] = gene

    def get(self, position: int) -> Gene:
        """Get a gene by its position

        :param position: Position of the gene in the contig

        :return: Wanted gene

        :raises TypeError: Position is not an integer
        """
        if not isinstance(position, int):
            raise TypeError(f"Position to get gene must be an integer. The provided type was {type(position)}")
        gene = self[position]
        if gene is None:
            logging.getLogger("PPanGGOLiN").debug("Given position result with a None Gene")
        return gene

    def remove(self, position):
        """Remove a gene by its position

        :param position: Position of the gene in the contig

        :raises TypeError: Position is not an integer
        """
        if not isinstance(position, int):
            raise TypeError(f"Position to get gene must be an integer. The provided type was {type(position)}")
        del self[position]

    def get_genes(self, begin: int = 0, end: int = None) -> List[Gene]:
        """
        Gets a list of genes within a range.
        If no arguments are given it return all genes.

        :param begin: Position of the first gene to retrieve
        :param end: Position of the last gene to not retrieve

        :return: List of genes between begin and end position

        :raises TypeError: If begin or end is not an integer
        :raises ValueError: If begin position is greater than end positon
        """

        if end is None:
            end = self.length

        if not isinstance(begin, int) or not isinstance(end, int):
            raise TypeError(f"Expected type int for 'begin' and 'end', but received types '{type(begin)}' and '{type(end)}'.")

        if begin >= end:
            raise ValueError("The 'begin' position must be less than the 'end' position.")

        else:
            return self._genes_position[begin: end]

    @property
    def number_of_genes(self) -> int:
        return len(self._genes_position)

    @property
    def genes(self) -> Generator[Gene, None, None]:
        """ Give the gene content of the contig

        :return: Generator of genes in contig
        """
        for gene in self._genes_position:
            if gene is not None:
                yield gene

    @property
    def organism(self) -> Organism:
        """Return organism that Feature belongs to.

        :return: Organism of the feature
        """
        return self._organism

    @organism.setter
    def organism(self, organism: Organism):
        """Set the organism belonging to the contig

        :param organism: Organism to set

        :raises TypeError: Given organism is not an instance Organism
        """
        if not isinstance(organism, Organism):
            raise TypeError(f'Expected type Organism, got {type(organism)}')
        self._organism = organism

    def add_rna(self, rna: RNA):
        """ Add RNA to contig

        :param rna: RNA object to add

        :raises TypeError: RNA is not instance RNA
        :raises KeyError: Another RNA with the same ID already exists in the contig
        """
        if not isinstance(rna, RNA):
            raise TypeError(f"'RNA' type was expected but you provided a '{type(rna)}' type object")
        if rna in self._rna_getter:
            raise KeyError(f"RNA with the id: {rna.ID} already exist in contig {self.name}")
        self._rna_getter.add(rna)

    @property
    def RNAs(self) -> Generator[RNA, None, None]:
        """Return all the RNA in the contig

        :return: Generator of RNA
        """
        yield from self._rna_getter

    @property
    def number_of_rnas(self) -> int:
        """Get the number of RNA in the contig
        """
        return len(self._rna_getter)


class Organism(MetaFeatures):
    """
    Describe the Genome content and some information

    Methods:
    - `families`: Returns a set of gene families present in the organism.
    - `genes`: Returns a generator to get genes in the organism.
    - `number_of_genes`: Returns the number of genes in the organism.
    - `contigs`: Returns the values in the contig dictionary from the organism.
    - `get_contig`: Gets the contig with the given identifier in the organism, adding it if it does not exist.
    - `_create_contig`: Creates a new contig object and adds it to the contig dictionary.
    - `mk_bitarray`: Produces a bitarray representing the presence/absence of gene families in the organism using the provided index.

    Fields:
    - `name`: Name of the organism.
    - `bitarray`: Bitarray representing the presence/absence of gene families in the organism.
    """

    def __init__(self, name: str):
        """Constructor Method

        :param name: Name of the genome
        """
        assert isinstance(name, str), "Organism name should be a string"
        assert name != "", "Organism name should not be empty"

        super().__init__()
        self.name = name
        self._contigs_getter = {}
        self._families = None
        self.bitarray = None

    def __str__(self):
        return self.name

    def _set_families(self):
        """Set the set of gene families belonging to organism
        """
        self._families = {gene.family for gene in self.genes}

    def __setitem__(self, name: str, contig: Contig):
        """ Set contig to the organism

        :param name: Name of the contig
        :param contig: Contig object to add in the organism

        :raises TypeError: If the contig is not instance Contig
        :raises TypeError: If the name is not instance string
        :raises KeyError: Contig with the given name already exist in the organism
        """

        if not isinstance(name, str):
            raise TypeError(f"Contig name should be a string. You provided a '{type(name)}' type object")
        if not isinstance(contig, Contig):
            raise TypeError(f"'Contig' type was expected but you provided a '{type(contig)}' type object")
        if name in self._contigs_getter:  # Add test if contig are equivalent when __eq__ method will be defined in Contig
            raise KeyError(f"Contig {contig.name} already in organism {self.name}")
        self._contigs_getter[contig.name] = contig
        contig.organism = self

    def __getitem__(self, name: str) -> Contig:
        """Get the contig for the given position

        :param name: Name of the contig

        :return:  Wanted contig for the given name

        :raises TypeError: If name is not a string
        :raises KeyError: Name does not exist in the organism
        """
        if not isinstance(name, str):
            raise TypeError(f"Expected type is string, given type was '{type(name)}'")
        try:
            return self._contigs_getter[name]
        except KeyError:
            raise KeyError(f"Contig with the name: {name} does not exist in the organism")

    def __delitem__(self, name):
        """Remove the contig for the given name

        :param name: Name of the contig

        :raises TypeError: If name is not a string
        :raises KeyError: Name does not exist in the organism
        """
        if not isinstance(name, int):
            raise TypeError(f"Expected type is int, given type was '{type(name)}'")
        try:
            del self._contigs_getter[name]
        except KeyError:
            raise KeyError("Position of the gene in the contig does not exist")

    def __len__(self):
        """ Get number of contigs in organism

        :return: Number of contigs in organism
        """
        return len(self._contigs_getter.keys())

    @property
    def families(self):
        """Return the gene families present in the organism

        :return: Generator of gene families
        :rtype: Generator[GeneFamily, None, None]
        """
        if self._families is None:
            self._set_families()
        yield from self._families

    def number_of_families(self) -> int:
        """Get the number of gene families in the organism

        :return: Number of gene families
        """
        if self._families is None:
            self._set_families()
        return len(self._families)

    @property
    def genes(self) -> Generator[Gene, None, None]:
        """Generator to get genes in the organism

        :return: Generator of genes
        """
        for contig in self.contigs:
            yield from contig.genes

    @property
    def rna_genes(self) -> Generator[RNA, None, None]:
        """Generator to get genes in the organism

        :return: Generator of genes
        """
        for contig in self.contigs:
            yield from contig.RNAs

    def number_of_genes(self) -> int:
        """ Get number of genes in the organism

        :return: Number of genes
        """
        return sum((contig.number_of_genes for contig in self.contigs))
    
    def number_of_rnas(self) -> int:
        """ Get number of genes in the organism

        :return: Number of genes
        """
        return sum((contig.number_of_rnas for contig in self.contigs))

    @property
    def contigs(self) -> Generator[Contig, None, None]:
        """ Generator of contigs in the organism

        :return: Values in contig dictionary from organism
        """
        yield from self._contigs_getter.values()


    @property
    def number_of_contigs(self) -> int:
        """ Get number of contigs in organism

        :return: Number of contigs in organism
        """
        return len(self._contigs_getter)


    def add(self, contig: Contig):
        """Add a contig to organism

        :param: Contig to add in organism

        :raises KeyError: Contig with the given name already exist in the organism
        """
        assert isinstance(contig, Contig), f"Contig object is expected, given type was {type(contig)}"
        try:
            _ = self.get(contig.name)
        except KeyError:
            self[contig.name] = contig
        else:
            raise KeyError(f"Contig {contig.name} already in organism {self.name}")

 
    def get(self, name: str) -> Contig:
        """
        Get contig with the given identifier in the organism

        :param name: Contig identifier

        :return: The contig with the given identifier
        """
        return self[name]


    def remove(self, name: str) -> Contig:
        """
        Remove a contig with the given identifier in the organism

        :param name: Contig identifier

        :return: The contig with the given identifier
        """
        del self[name]


    def mk_bitarray(self, index: Dict[Organism, int], partition: str = 'all'):
        """Produces a bitarray representing the presence / absence of families in the organism using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.

        :param partition: Filters partition
        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`

        :raises Exception: Partition is not recognized
        """
        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == 'all':
            logging.getLogger("PPanGGOLiN").debug("all")
            for fam in self.families:
                self.bitarray[index[fam]] = 1
        elif partition in ['shell', 'cloud']:
            logging.getLogger("PPanGGOLiN").debug("shell, cloud")
            for fam in self.families:
                if fam.named_partition == partition:
                    self.bitarray[index[fam]] = 1
        elif partition == 'accessory':
            logging.getLogger("PPanGGOLiN").debug("accessory")
            for fam in self.families:
                if fam.named_partition in ['shell', 'cloud']:
                    self.bitarray[index[fam]] = 1
        else:
            raise Exception("There is not any partition corresponding please report a github issue")
        
