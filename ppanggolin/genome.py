#!/usr/bin/env python3

from __future__ import annotations
import logging
from typing import Dict, Generator, List, Union, Set, Tuple, Iterable
from collections import defaultdict

# installed libraries
import gmpy2

# local libraries
from ppanggolin.metadata import MetaFeatures
from ppanggolin.utils import get_consecutive_region_positions


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
        if identifier == "":
            raise ValueError("Identifier should not be empty")
        super().__init__()
        self.ID = identifier
        self.is_fragment = False
        self.type = ""
        self.start = None
        self.stop = None
        self.coordinates = None
        self.strand = None
        self.product = None
        self.name = None
        self.local_identifier = None
        self._organism = None
        self._contig = None
        self.dna = None

    def __str__(self) -> str:
        """String representation of the feature

        :return: feature identifier
        """
        return str(self.ID)

    def __len__(self) -> int:
        """Return gene length

        :return: gene length

        :raises ValueError: If coordinates are not defined in gene
        """

        try:
            return sum([(stop - start + 1) for start, stop in self.coordinates])
        except TypeError:
            raise ValueError(
                f"Coordinates of gene {self} have not been defined. Getting its length is then impossible."
            )

    @property
    def has_joined_coordinates(self) -> bool:
        """
        Whether or not the feature has joined coordinates.

        """
        if len(self.coordinates) > 1:
            return True
        else:
            return False

    @property
    def overlaps_contig_edge(self) -> bool:
        """
        Check based on the coordinates of the feature, if the gene seems to overlap contig edge.

        """

        start_stop = self.coordinates[0]
        for start_stop_next in self.coordinates[1:]:
            if start_stop > start_stop_next:
                return True
            start_stop = start_stop_next

        return False

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
            raise TypeError(f"Expected type Organism, got {type(organism)}")
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
            raise TypeError(f"Expected type Contig, got {type(contig)}")
        self._contig = contig

    def fill_annotations(
        self,
        start: int,
        stop: int,
        strand: str,
        gene_type: str = "",
        name: str = "",
        product: str = "",
        local_identifier: str = "",
        coordinates: List[Tuple[int, int]] = None,
    ):
        """
        Fill general annotation for child classes

        :param start: Start position
        :param stop: Stop position
        :param coordinates: start and stop positions. in a list of tuple. Can have multiple tuple in case of join gene
        :param strand: associated strand
        :param gene_type: Type of gene
        :param name: Name of the feature
        :param product: Associated product
        :param local_identifier: Identifier provided by the original file

        :raises TypeError: If attribute value does not correspond to the expected type
        :raises ValueError: If strand is not '+' or '-'
        """
        if coordinates is None:
            coordinates = [(start, stop)]

        if not isinstance(start, int):
            raise TypeError(
                f"Start should be int. Got {type(start)} instead in {self} from {self.organism}."
            )
        if not isinstance(stop, int):
            raise TypeError(
                f"Stop should be int. Got {type(stop)} instead in {self} from {self.organism}."
            )
        if not isinstance(strand, str):
            raise TypeError(
                f"Strand should be str. Got {type(strand)} instead in {self} from {self.organism}."
            )
        if not isinstance(gene_type, str):
            raise TypeError(
                f"Gene type should be str. Got {type(gene_type)} instead in {self} from {self.organism}."
            )
        if not isinstance(name, str):
            raise TypeError(
                f"Name should be str. Got {type(name)} instead in {self} from {self.organism}."
            )
        if not isinstance(product, str):
            raise TypeError(
                f"Product should be str. Got {type(product)} instead in {self} from {self.organism}."
            )
        if not isinstance(local_identifier, str):
            raise TypeError(
                f"Local identifier should be str. Got {type(local_identifier)} instead in {self} from {self.organism}."
            )
        if strand not in ["+", "-"]:
            raise ValueError(
                f"Strand should be '+' or '-'. Got {strand} instead in {self} from {self.organism}."
            )
        if not isinstance(coordinates, list):
            raise TypeError(
                f"Coordinates should be of type list. Got {type(coordinates)} instead in {self} from {self.organism}."
            )

        for start_i, stop_i in coordinates:
            if not isinstance(start_i, int):
                raise TypeError(
                    f"Start should be int. Got {type(start_i)} instead in {self} from {self.organism}."
                )
            if not isinstance(stop_i, int):
                raise TypeError(
                    f"Stop should be int. Got {type(stop_i)} instead in {self} from {self.organism}."
                )
            if stop_i < start_i:
                raise ValueError(
                    f"Wrong coordinates: {coordinates}. Start ({start_i}) should not be greater than stop ({stop_i}) in {self} from {self.organism}."
                )
            if start_i < 1 or stop_i < 1:
                raise ValueError(
                    f"Wrong coordinates: {coordinates}. Start ({start_i}) and stop ({stop_i}) should be greater than 0 in {self} from {self.organism}."
                )

        self.start = start
        self.stop = stop
        self.strand = strand
        self.type = gene_type
        self.product = product
        self.name = name
        self.local_identifier = local_identifier
        self.coordinates = coordinates

    def fill_parents(self, organism: Organism = None, contig: Contig = None):
        """Associate object to an organism and a contig

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
        assert isinstance(
            sequence, str
        ), f"'str' type was expected for dna sequence but you provided a '{type(sequence)}' type object"

        self.dna = sequence

    def string_coordinates(self) -> str:
        """
        Return a string representation of the coordinates
        """
        return ",".join(f"{start}..{stop}" for start, stop in self.coordinates)

    def start_relative_to(self, gene):
        """ """
        if gene.start <= self.start:
            return self.start
        if gene.start > self.start:
            return self.start + self.contig.length

    def stop_relative_to(self, gene):
        """ """
        if gene.start <= self.stop:
            return self.stop

        if gene.start > self.stop:
            return self.stop + self.contig.length


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
        self.is_partial = False  # is the gene a partial gene ?
        self._frame = None  # One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..

    @property
    def family(self):
        """Return GeneFamily that Gene belongs to.

        :return: Gene family of the gene
        :rtype: GeneFamily
        """
        return self._family

    @family.setter
    def family(self, family):
        """Set the GeneFamily belonging to the gene

        :param family: Gene family linked to the gene
        """
        from ppanggolin.geneFamily import GeneFamily

        if not isinstance(family, GeneFamily):
            raise TypeError(f"Expected type GeneFamily, got {type(family)}")
        self._family = family

    @property
    def RGP(self):
        """Return the RGP that gene belongs to

        :return: RGP of the Gene
        :rtype: Region
        """
        return self._RGP

    @RGP.setter
    def RGP(self, region):
        """Set the Region belonging to the gene

        :param region: Region linked to the gene
        """
        from ppanggolin.region import Region

        if not isinstance(region, Region):
            raise TypeError(f"Expected type Organism, got {type(region)}")
        self._RGP = region

    @property
    def spot(self):
        """Get the spot belonging to the gene

        :return: the spot linked to the gene
        :rtype: Spot
        """
        if self.RGP is not None and self.RGP.spot is not None:
            return self.RGP.spot
        else:
            return None

    @property
    def module(self):
        """
        Get the modules belonging to the gene

        :return: get the modules linked to the gene
        :rtype: Module
        """
        return self.family.module

    def fill_annotations(
        self,
        position: int = None,
        genetic_code: int = 11,
        is_partial: bool = False,
        frame: int = 0,
        **kwargs,
    ):
        """Fill Gene annotation provide by PPanGGOLiN dependencies

        :param position: Gene localization in genome
        :param genetic_code: Genetic code associated to gene
        :param is_partial: is the gene a partial gene
        :param frame: One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon,
                      '1' that the second base is the first base of a codon, and so on..
        :param kwargs: look at Feature.fill_annotations methods

        :raises TypeError: If position or genetic code value is not instance integers
        """
        super().fill_annotations(**kwargs)
        if position is not None and not isinstance(position, int):
            raise TypeError("position should be an integer")
        if not isinstance(genetic_code, int):
            raise TypeError("Genetic code should be an integer")

        if not isinstance(is_partial, bool):
            raise TypeError("partial code should be an boolean")

        self.position = position
        self.genetic_code = genetic_code
        self.is_partial = is_partial
        self.frame = frame

    def add_protein(self, protein: str):
        """Add a protein sequence corresponding to translated gene

        :param protein: Protein sequence

        :raise TypeError: Protein sequence must be a string
        """
        if not isinstance(protein, str):
            raise TypeError(
                f"'str' type was expected but you provided a '{type(protein)}' type object"
            )
        self.protein = protein

    @property
    def frame(self) -> int:
        """
        Get the frame of the gene

        """
        assert (
            self._frame is not None
        ), "frame is already set and should not be set another time."

        return self._frame

    @frame.setter
    def frame(self, frame: int):
        """Set the length of the contig

        :param contig_len: length of the contig
        """
        assert (
            self._frame is None
        ), "frame is already set and should not be set another time."

        if frame not in [0, 1, 2]:
            raise ValueError("Frame should be equal to 0, 1 or 2.")

        self._frame = frame


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

    TODO: Getter gene should be based on gene ID, and 2 other attributes should exist to get them by start or position.
          Also, when set a new gene in contig, start, stop and strand should be check to check difference, maybe define __eq__ method in gene class.
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
        self._rna_getter = (
            set()
        )  # Saving the rna annotations. We're not using them in the vast majority of cases.
        self._genes_getter = {}
        self._genes_position = []
        self._organism = None
        self._length = None

    def __str__(self) -> str:
        """Returns a string representation of the contig

        :return: Name of the contig
        """
        return self.name

    def __setitem__(self, coordinate: Tuple[int, int, str], gene: Gene):
        """
        Set gene to Contig

        Check if a gene with the same coordinate exists already in the contig.

        :param coordinate: Tuple containing start, stop and strand of the gene
        :param gene: Gene object to add

        :raises TypeError: If the gene is not instance Gene
        :raises ValueError: If a gene in getter already exists at the start
        :raises AttributeError: If the gene position in the contig is not fill
        """

        if not isinstance(gene, Gene):
            raise TypeError(
                f"'Gene' type was expected but you provided a '{type(gene)}' type object"
            )

        if coordinate in self._genes_getter:
            raise ValueError(
                f"Gene '{self._genes_getter[coordinate].ID}' with coordinate {coordinate} already exists in the "
                f"contig '{self.name}' {f'from genome {self.organism}' if self.organism else ''}, "
                f"cannot add gene '{gene.ID}' {f'from genome {gene.organism}' if gene.organism else ''}"
            )

        if gene.position is None:
            raise AttributeError(
                "The gene object needs to have its position in the contig filled before adding it"
            )

        # Adding empty values.
        # They should be filled by the end of the parsing.
        # Doing this because genes are not always met in order.
        self._genes_position.extend(
            [None] * (gene.position - len(self._genes_position) + 1)
        )
        self._genes_position[gene.position] = gene
        self._genes_getter[coordinate] = gene

    # TODO define eq function

    @property
    def length(self) -> Union[int, None]:
        """Get the length of the contig"""
        return self._length

    @length.setter
    def length(self, contig_len: int):
        """Set the length of the contig

        :param contig_len: length of the contig
        """
        if not isinstance(contig_len, int):
            raise TypeError("Contig length is expected to be an integer")
        if contig_len < 0:
            raise ValueError("Contig length must be positive")

        if self._length is None:
            self._length = contig_len
        elif self.length != contig_len:
            logging.getLogger("PPanGGOLiN").debug(
                f"Known contig length = {self.length}, new length = {contig_len}"
            )
            raise ValueError(
                "Attempting to define a contig length different from the previously defined value."
            )

    def __len__(self) -> int:
        """Get the length of the contig

        :return: contig length

        :todo: It could be better that len return the number of genes in the contig
        """
        if self.length is None:
            raise ValueError(
                f"Contig length of {self.name} from genome {self.organism} has not been defined. Getting its length is then impossible."
            )
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
            gene = self._genes_position[position]
        except KeyError:
            raise KeyError("Position of the gene in the contig does not exist")
        else:
            return gene

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
            raise TypeError(
                f"Unexpected class / type for {type(gene)} when adding it to a contig"
            )

        for attr in ["start", "stop", "position", "strand"]:
            if getattr(gene, attr) is None:
                raise AttributeError(f"Gene {gene.name} is not fill with {attr}")

        if gene.strand not in ["+", "-"]:
            raise AttributeError(
                f"Strand of Gene {gene.name} does not have the expected format. Expect '-' or '+' got {gene.strand}"
            )

        self[(gene.start, gene.stop, gene.strand)] = gene

    def get_by_coordinate(self, coordinate: Tuple[int, int, str]) -> Gene:
        """
        Get a gene by its coordinate

        :param coordinate: Tuple containing start, stop and strand of the gene

        :return: The gene with the specified coordinate.

        :raises TypeError: Position is not an integer
        """
        if not isinstance(coordinate, Tuple):
            raise TypeError(
                f"Coordinate to get gene must be a tuple. The provided type was {type(coordinate)}"
            )

        gene = self[coordinate]
        if gene is None:
            logging.getLogger("PPanGGOLiN").debug(
                "Given position result with a None Gene"
            )
        return gene

    def remove(self, position):
        """Remove a gene by its position

        :param position: Position of the gene in the contig

        :raises TypeError: Position is not an integer
        """
        if not isinstance(position, int):
            raise TypeError(
                f"Position to get gene must be an integer. The provided type was {type(position)}"
            )
        del self[position]

    def get_genes(
        self, begin: int = 0, end: int = None, outrange_ok: bool = False
    ) -> List[Gene]:
        """
        Gets a list of genes within a range of gene position.
        If no arguments are given it return all genes.

        :param begin: Position of the first gene to retrieve
        :param end: Position of the last gene to not retrieve
        :param outrange_ok: If True even is the last position is out of range return all the genes from begin to last position

        :return: List of genes between begin and end position

        :raises TypeError: If begin or end is not an integer
        :raises ValueError: If begin position is greater than end position
        :raises IndexError: If end position is greater than last gene position in contig
        """
        if end is None:
            end = self._genes_position[-1].position

        if not isinstance(begin, int) or not isinstance(end, int):
            raise TypeError(
                f"Expected type int for 'begin' and 'end', "
                f"but received types '{type(begin)}' and '{type(end)}'."
            )

        if begin > end:
            raise ValueError(
                "The 'begin' position must be less than the 'end' position."
            )

        if end > self._genes_position[-1].position:
            if outrange_ok:
                end = self._genes_position[-1].position
            else:
                raise IndexError(f"Gene at position {end} is out of range")

        if end == self._genes_position[-1].position:
            return self._genes_position[begin:]
        else:
            if begin == end:
                return self._genes_position[begin]
            else:
                return self._genes_position[begin:end]

    @property
    def number_of_genes(self) -> int:
        """Get the number of genes in the contig

        :return: the number of genes in the contig
        """
        return len(self._genes_position)

    @property
    def genes(self) -> Generator[Gene, None, None]:
        """Give the gene content of the contig

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
            raise TypeError(f"Expected type Organism, got {type(organism)}")
        self._organism = organism

    def add_rna(self, rna: RNA):
        """Add RNA to contig

        :param rna: RNA object to add

        :raises TypeError: RNA is not instance RNA
        :raises KeyError: Another RNA with the same ID already exists in the contig
        """
        if not isinstance(rna, RNA):
            raise TypeError(
                f"'RNA' type was expected but you provided a '{type(rna)}' type object"
            )
        if rna in self._rna_getter:
            raise KeyError(
                f"RNA with the id: {rna.ID} already exist in contig {self.name}"
            )
        self._rna_getter.add(rna)

    @property
    def RNAs(self) -> Generator[RNA, None, None]:
        """Return all the RNA in the contig

        :return: Generator of RNA
        """
        yield from self._rna_getter

    @property
    def number_of_rnas(self) -> int:
        """Get the number of RNA in the contig"""
        return len(self._rna_getter)

    def add_contig_length(self, contig_length: int):
        """
        Add contig length to Contig object.

        :param contig_length: Length of the contig.
        :raises ValueError: If trying to define a contig length different from previously defined.
        """
        if self.length is None:
            self.length = contig_length

        elif self.length != contig_length:
            raise ValueError(
                "Attempting to define a contig length different from the previously defined value."
            )

    @property
    def regions(self):
        """Get the regions belonging to this contig

        :return: RGP in the contig
        :rtype: Generator[Region, None, None]
        """
        regions = set()
        for gene in self.genes:
            if gene.RGP is not None:
                regions.add(gene.RGP)
        yield from regions

    @property
    def spots(self):
        """Get the spot belonging to this contig

        :return: Spot in the contig
        :rtype: Generator[Spot, None, None]
        """
        spots = set()
        for region in self.regions:
            if region.spot is not None:
                spots.add(region.spot)
        yield from spots

    @property
    def families(self):
        """Get the families belonging to this contig

        :return: families in the contig
        :rtype: Generator[GeneFamily, None, None]
        """
        families = set()
        for gene in self.genes:
            if gene.family is None:
                raise ValueError(
                    "Gene has no family, that should not happen. "
                    "Check if you're families has been computed or loaded."
                    "If it's the case, you can report an issue on our GitHub."
                )
            families.add(gene.family)
        yield from families

    @property
    def modules(self):
        """Get the modules belonging to this contig

        :return: Modules belonging to this contig
        :rtype: Generator[Module, None, None]
        """
        modules = set()
        for family in self.families:
            for module in family.modules:
                modules.add(module)
        yield from modules

    def get_ordered_consecutive_genes(self, genes: Iterable[Gene]) -> List[List[Gene]]:
        """
        Order the given genes considering the circularity of the contig.

        :param genes: An iterable containing genes supposed to be consecutive along the contig.
        :return: A list of lists containing ordered consecutive genes considering circularity.
        """
        gene_positions = [gene.position for gene in genes]

        # Determine consecutive region positions
        consecutive_region_positions = get_consecutive_region_positions(
            region_positions=gene_positions, contig_gene_count=self.number_of_genes
        )

        consecutive_genes_lists = [
            [self[position] for position in consecutive_positions]
            for consecutive_positions in consecutive_region_positions
        ]

        return consecutive_genes_lists


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
    - `mk_bitarray`: Produces a bitarray representing the presence/absence of gene families
      in the organism using the provided index.

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

    def __str__(self) -> str:
        """String representation of the genome

        :return: Name of the genome
        """
        return self.name

    def _set_families(self):
        """Set the set of gene families belonging to organism"""
        self._families = {gene.family for gene in self.genes}

    def __setitem__(self, name: str, contig: Contig):
        """Set contig to the organism

        :param name: Name of the contig
        :param contig: Contig object to add in the organism

        :raises TypeError: If the contig is not instance Contig
        :raises TypeError: If the name is not instance string
        :raises KeyError: Contig with the given name already exist in the organism
        """

        if not isinstance(name, str):
            raise TypeError(
                f"Contig name should be a string. You provided a '{type(name)}' type object"
            )
        if not isinstance(contig, Contig):
            raise TypeError(
                f"'Contig' type was expected but you provided a '{type(contig)}' type object"
            )
        if name in self._contigs_getter:
            # Add test if contig are equivalent when __eq__ method will be defined in Contig
            raise KeyError(f"Contig {contig.name} already in genome {self.name}")
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
            contig = self._contigs_getter[name]
        except KeyError:
            raise KeyError(f"Contig with the name: {name} does not exist in the genome")
        else:
            return contig

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
        """Get number of contigs in organism

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
        """Get number of genes in the organism

        :return: Number of genes
        """
        return sum(contig.number_of_genes for contig in self.contigs)

    def number_of_rnas(self) -> int:
        """Get number of genes in the organism

        :return: Number of genes
        """
        return sum(contig.number_of_rnas for contig in self.contigs)

    @property
    def contigs(self) -> Generator[Contig, None, None]:
        """Generator of contigs in the organism

        :return: Values in contig dictionary from organism
        """
        yield from self._contigs_getter.values()

    @property
    def number_of_contigs(self) -> int:
        """Get number of contigs in organism

        :return: Number of contigs in organism
        """
        return len(self._contigs_getter)

    def add(self, contig: Contig):
        """Add a contig to organism

        :param: Contig to add in organism

        :raises KeyError: Contig with the given name already exist in the organism
        """
        assert isinstance(
            contig, Contig
        ), f"Contig object is expected, given type was {type(contig)}"
        try:
            _ = self.get(contig.name)
        except KeyError:
            self[contig.name] = contig
        else:
            raise KeyError(f"Contig {contig.name} already in genome {self.name}")

    def get(self, name: str) -> Contig:
        """
        Get contig with the given identifier in the organism

        :param name: Contig identifier

        :return: The contig with the given identifier
        """
        return self[name]

    def remove(self, name: str):
        """
        Remove a contig with the given identifier in the organism

        :param name: Contig identifier
        """
        del self[name]

    @property
    def modules(self):
        """
        Get all the modules belonging to this genome

        :return: Generator of modules
        :rtype: Generator[Module, None, None]
        """
        modules = {family.module for family in self.families if family.has_module}
        yield from modules

    @property
    def number_of_modules(self) -> int:
        """
        Get number of modules in organism

        :return: Number of modules in organism
        """
        return len(list(self.modules))

    @property
    def regions(self):
        """Get all RGPS belonging to this genome

        :return: Generator of RGPS
        :rtype: Generator[Region, None, None]
        """
        regions = set()
        for gene in self.genes:
            if gene.RGP is not None:
                regions.add(gene.RGP)
        yield from regions

    @property
    def number_of_regions(self) -> int:
        """
        Get number of RGP in organism

        :return: Number of RGP in organism
        """
        return len(list(self.regions))

    @property
    def spots(self):
        """Get all spots belonging to this genome

        :return: Generator of spots
        :rtype: Generator[Spot, None, None]
        """
        spots = set()
        for region in self.regions:
            if region.spot is not None:
                spots.add(region.spot)
        yield from spots

    @property
    def number_of_spots(self) -> int:
        """
        Get number of spots in organism

        :return: Number of spots in organism
        """
        return len(list(self.spots))

    def mk_bitarray(self, index: Dict[Organism, int], partition: str = "all"):
        """Produces a bitarray representing the presence / absence of families in the organism using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.

        :param partition: Filters partition
        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`

        :raises Exception: Partition is not recognized
        """
        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == "all":
            logging.getLogger("PPanGGOLiN").debug("all")
            for fam in self.families:
                self.bitarray[index[fam]] = 1
        elif partition in ["shell", "cloud"]:
            logging.getLogger("PPanGGOLiN").debug("shell, cloud")
            for fam in self.families:
                if fam.named_partition == partition:
                    self.bitarray[index[fam]] = 1
        elif partition == "accessory":
            logging.getLogger("PPanGGOLiN").debug("accessory")
            for fam in self.families:
                if fam.named_partition in ["shell", "cloud"]:
                    self.bitarray[index[fam]] = 1
        else:
            raise ValueError(
                "There is not any partition corresponding please report a github issue"
            )

    def group_genes_by_partition(self) -> Dict[str, Set]:
        """
        Groups genes based on their family's named partition and returns a dictionary
        mapping partition names to sets of genes belonging to each partition.

        :return: A dictionary containing sets of genes grouped by their family's named partition.
        """
        partition_to_gene = defaultdict(set)
        contigs_count = 0

        for contig in self.contigs:
            contigs_count += 1
            for gene in contig.genes:
                partition_to_gene[gene.family.named_partition].add(gene)

        return partition_to_gene
