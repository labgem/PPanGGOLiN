#!/usr/bin/env python3

# default libraries
from __future__ import annotations
import logging

# installed libraries
import networkx as nx
from typing import Dict, Generator, List, Set, Union, Tuple
import gmpy2

# local libraries
from ppanggolin.genome import Gene, Organism, Contig
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.metadata import MetaFeatures
from ppanggolin.utils import (
    find_region_border_position,
    get_consecutive_region_positions,
)


class Region(MetaFeatures):
    """
    The 'Region' class represents a region of genomic plasticity.

    Methods:
        - 'genes': the property that generates the genes in the region as they are ordered in contigs.
        - 'families': the property that generates the gene families in the region.
        - 'Length': the property that gets the length of the region.
        - 'organism': the property that gets the organism linked to the region.
        - 'Contig': the property that gets the starter contig linked to the region.
        - 'is_whole_contig': the property that indicates if the region is an entire contig.
        - 'is_contig_border': the property that indicates if the region is bordering a contig.
        - 'get_rnas': the method that gets the RNA in the region.
        - 'Get_bordering_genes': the method that gets the bordered genes in the region.

    Fields:
        - 'name': the name of the region.
        - 'score': the score of the region.
        - 'Starter': the first gene in the region.
        - 'stopper': the last gene in the region.
    """

    id_counter = 0

    def __init__(self, name: str):
        """Constructor method

        :param name: Name of the region
        """
        super().__init__()
        self._genes_getter = {}
        self.name = name
        self.score = None
        self._starter = None
        self._stopper = None
        self._coordinates = None
        self._overlaps_contig_edge = None
        self._contig = None
        self._organism = None
        self.ID = Region.id_counter
        self._spot = None
        self.projected = False  # If the rgp is from a projected genome. If true can have multiple spots
        Region.id_counter += 1

    def __str__(self):
        return self.name

    def __repr__(self) -> str:
        """
        Region representation
        """
        return f"RGP name:{self.name}"

    def __hash__(self) -> int:
        """Create a hash value for the region"""
        return id(self)

    def __lt__(self, obj):
        return self.ID < obj.ID

    def __gt__(self, obj):
        return self.ID > obj.ID

    def __eq__(self, other: Region) -> bool:
        """
        Test whether two Region objects have the same gene families

        :param other: Another region to test equality of regions

        :return: Equal or not

        :raises TypeError: Try to compare a region with another type object
        """
        if not isinstance(other, Region):
            raise TypeError(
                f"'Region' type object was expected, but '{type(other)}' type object was provided."
            )
        if [gene.family for gene in self.genes] == [
            gene.family for gene in other.genes
        ]:
            return True
        if [gene.family for gene in self.genes] == [
            gene.family for gene in list(other.genes)[::-1]
        ]:
            return True
        return False

    def __len__(self) -> int:
        """Get the number of genes in the region"""
        return len(self._genes_getter)

    def __setitem__(self, position: int, gene: Gene):
        """Set a gene by is position in the region

        :param position: Position of the gene in the contig
        :param gene: Gene to add in the region

        :raises TypeError: If the gene is not an instance of Gene.
        :raises ValueError: If the organism or contig of the gene is different from the region.
        :raises KeyError: If another gene already exists at the specified position.
        :raises ValueError: If the position of the gene does not match the provided position.
        """

        if position != gene.position:
            raise ValueError(
                f"The given gene position ({position}) to set the gene in the region and the position of the gene ({gene.position})  are different. "
            )

        if len(self) == 0:
            # first gene to be added to the region
            self._organism = gene.organism
            self._contig = gene.contig

        if len(self) > 0:
            if gene.organism != self.organism:
                raise ValueError(
                    f"Gene {gene.name} is from a different genome than the first defined in RGP. "
                    "That's not possible"
                )
            if gene.contig != self.contig:
                raise ValueError(
                    f"Gene {gene.name} is from a different contig than the first defined in RGP. "
                    "That's not possible"
                )
        if position in self._genes_getter and self[position] != gene:
            raise KeyError("Another gene already exist at this position")
        self._genes_getter[position] = gene

        # Adding a new gene imply to reidentify first (starter) and last (stopper) genes of the rgp.
        self._starter = None
        self._stopper = None
        self._coordinates = None
        self._overlaps_contig_edge = None

        gene.RGP = self

    def identify_rgp_last_and_first_genes(self):
        """
        Identify first and last genes of the rgp by taking into account the circularity of contigs.

        Set the attributes _starter: first gene of the region  and _stopper: last gene of the region and _coordinates

        """
        rgp_genes_positions = list(self._genes_getter.keys())

        if len(rgp_genes_positions) == 0:
            raise ValueError(f"RGP ({self.name}) has no gene associated.")

        gene = self._genes_getter[rgp_genes_positions[0]]  # get a gene of the region
        first_gene_position, last_gene_position = find_region_border_position(
            region_positions=rgp_genes_positions,
            contig_gene_count=gene.contig.number_of_genes,
        )

        self._starter = self._genes_getter[first_gene_position]
        self._stopper = self._genes_getter[last_gene_position]

        if self._starter.start > self._stopper.stop:
            # this means region is overlapping the contig edge
            if not gene.contig.is_circular:
                raise ValueError(
                    f"Region seems to be overlapping the contig (first gene {self._starter.position}:{self._starter.coordinates} "
                    f"and last gene {self._stopper.position}:{self._stopper.coordinates} ) "
                    f"but the contig is not circular. This is unexpected. {rgp_genes_positions}"
                )

            self._coordinates = [
                (self._starter.start, self._starter.contig.length),
                (1, self._stopper.stop),
            ]
            self._overlaps_contig_edge = True
        else:
            self._coordinates = [(self._starter.start, self._stopper.stop)]
            self._overlaps_contig_edge = False

    def get_ordered_genes(self) -> List[Gene]:
        """
        Get ordered genes of the region, taking into account the circularity of contigs.

        :return: A list of genes ordered by their positions in the region.
        """

        rgp_genes_positions = list(self._genes_getter.keys())

        gene = self._genes_getter[rgp_genes_positions[0]]  # get a gene of the region

        consecutive_region_positions = get_consecutive_region_positions(
            region_positions=rgp_genes_positions,
            contig_gene_count=gene.contig.number_of_genes,
        )

        ordered_genes = [
            self._genes_getter[position]
            for ordered_positions in consecutive_region_positions
            for position in ordered_positions
        ]

        return ordered_genes

    def __getitem__(self, position: int) -> Gene:
        """Get the gene at the given position

        :param position: Position of the gene

        :return: Gene in the Region at the given position

        :raises KeyError: Gene at the given position does not exist
        """
        try:
            gene = self._genes_getter[position]
        except KeyError:
            raise KeyError(
                f"There is no gene at position {position} in RGP {self.name}"
            )
        else:
            return gene

    @property
    def starter(self) -> Gene:
        """
        Return first gene of the region. If this gene is not identified, it does that first.
        :return:  first gene of the region
        """
        if self._starter is None:
            self.identify_rgp_last_and_first_genes()

        return self._starter

    @property
    def stopper(self) -> Gene:
        """
        Return last gene of the region. If this gene is not identified, it does that first.
        :return: last gene of the region
        """
        if self._stopper is None:
            self.identify_rgp_last_and_first_genes()
        return self._stopper

    @property
    def coordinates(self) -> List[Tuple[int]]:
        """
        Return the coordinates of the region
        :return: coordinates of the region
        """
        if self._coordinates is None:
            self.identify_rgp_last_and_first_genes()
        return self._coordinates

    def string_coordinates(self) -> str:
        """
        Return a string representation of the coordinates
        """
        return ",".join([f"{start}..{stop}" for start, stop in self.coordinates])

    @property
    def overlaps_contig_edge(self) -> bool:
        if self._overlaps_contig_edge is None:
            self.identify_rgp_last_and_first_genes()
        return self._overlaps_contig_edge

    @property
    def spot(self) -> Union[Spot, None]:
        return self._spot

    @spot.setter
    def spot(self, spot: Spot):
        """Sets the spot of the RGP

        :param spot: spot to which the RGP is added

        :raise TypeError: if the given spot is not a Spot.
        """
        if isinstance(spot, Spot):
            self._spot = spot  # only 1 spot possible
        else:
            raise TypeError(
                f"Unexpected class / type for {type(spot)} when adding spot to a RGP"
            )

    def __delitem__(self, position):
        """Remove the gene at the given position

        :param position: Position of the gene

        :raises KeyError: Gene at the given position does not exist"""
        try:
            del self._genes_getter[position]
        except KeyError:
            raise KeyError(
                f"There is no gene at position {position} in RGP {self.name}"
            )

    def add(self, gene: Gene):
        """Add a gene to the region

        :param gene: Gene to add
        """
        if not isinstance(gene, Gene):
            raise TypeError(
                f"Unexpected class / type for {type(gene)} "
                f"when adding it to a region of genomic plasticity"
            )
        if gene.position is None:
            raise AttributeError(f"Gene {gene.name} is not fill with position")
        self[gene.position] = gene

    def get(self, position: int) -> Gene:
        """Get a gene by its position

        :param position: Position of the gene in the contig

        :return: Wanted gene

        :raises TypeError: Position is not an integer
        """
        if not isinstance(position, int):
            raise TypeError(
                f"Position to get gene must be an integer. The provided type was {type(position)}"
            )
        return self[position]

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

    @property
    def genes(self) -> Generator[Gene, None, None]:
        """Generate the gene as they are ordered in contigs

        :return: Genes in the region
        """
        yield from sorted(self._genes_getter.values(), key=lambda x: x.position)

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """Get the gene families in the RGP

        :return: Gene families
        """
        for gene in self.genes:
            yield gene.family

    @property
    def modules(self) -> Set[Module]:
        """Get the modules of gene families in the RGP

        :return: Modules found in families of the RGP
        """
        modules = {
            family.module for family in self.families if family.module is not None
        }
        return modules

    @property
    def number_of_families(self) -> int:
        """Get the number of different gene families in the region

        :return: Number of families
        """
        return len(set(self.families))

    @property
    def length(self):
        """Get the length of the region

        :return: Size of the region
        """
        return sum([(stop - start + 1) for start, stop in self.coordinates])

    @property
    def organism(self) -> Organism:
        """Get the Organism link to RGP

        :return: Organism corresponding to the region
        """
        return self._organism

    @property
    def contig(self) -> Contig:
        """Get the starter contig link to RGP

        :return: Contig corresponding to the region
        """
        return self._contig

    @property
    def start(self) -> int:
        """
        Get the starter start link to RGP

        :return: start position in the contig of the first gene of the RGP
        """
        return self.starter.start

    @property
    def stop(self) -> int:
        """
        Get the stopper stop link to RGP

        :return: start position in the contig of the last gene of the RGP
        """
        return self.stopper.stop

    @property
    def is_whole_contig(self) -> bool:
        """Indicates if the region is an entire contig

        :return: True if whole contig else False
        """
        if (
            self.starter.position == 0
            and self.stopper.position == self.contig.number_of_genes - 1
        ):
            return True
        return False

    @property
    def is_contig_border(self) -> bool:
        """Indicates if the region is bordering a contig

        :return: True if bordering else False

        :raises AssertionError: No genes in the regions, it's not expected
        """
        assert len(self) > 0, "Your region has no genes. Something wrong happened."

        if not self.contig.is_circular:
            first_gene = self.contig[0]
            last_gene = self.contig[-1]

            if self.starter == first_gene or self.stopper == last_gene:
                return True
        return False

    def get_bordering_genes(
        self, n: int, multigenics: Set[GeneFamily], return_only_persistents: bool = True
    ) -> List[List[Gene], List[Gene]]:
        """
        Get the bordered genes in the region. Find the n persistent and single copy gene bordering the region.
        If return_only_persistents is False, the method return all genes included between the n single copy and persistent genes.

        :param n: Number of genes to get
        :param multigenics: pangenome graph multigenic persistent families
        :param return_only_persistents: return only non multgenic persistent genes identify as the region.
                                        If False return all genes included between
                                        the borders made of n persistent and single copy genes around the region.

        :return: A list of bordering genes in start and stop position
        """
        genes_in_region = list(self.genes)
        # Identifying left border
        left_border = []
        pos = self.starter.position
        init = pos
        single_copy_persistent_count = 0
        while single_copy_persistent_count < n and (
            pos != 0 or self.contig.is_circular
        ):
            curr_gene = None
            if pos == 0:
                if self.contig.is_circular:
                    curr_gene = self.contig[pos - 1]
            else:
                curr_gene = self.contig[pos - 1]

            if (
                curr_gene is not None
                and curr_gene.family not in multigenics
                and curr_gene.family.named_partition == "persistent"
                and curr_gene not in genes_in_region
            ):
                left_border.append(curr_gene)
                single_copy_persistent_count += 1
            elif (
                curr_gene is not None
                and curr_gene not in genes_in_region
                and not return_only_persistents
            ):
                left_border.append(curr_gene)

            pos -= 1
            if pos == -1 and self.contig.is_circular:
                pos = self.contig.number_of_genes
            if pos == init:
                break  # looped around the contig

        # Identifying right border
        right_border = []
        pos = self.stopper.position
        init = pos
        single_copy_persistent_count = 0
        while single_copy_persistent_count < n and (
            pos != self.contig.number_of_genes - 1 or self.contig.is_circular
        ):
            curr_gene = None
            if pos == self.contig.number_of_genes - 1:
                if self.contig.is_circular:
                    curr_gene = self.contig[0]
            else:
                curr_gene = self.contig[pos + 1]
            if (
                curr_gene is not None
                and curr_gene.family not in multigenics
                and curr_gene.family.named_partition == "persistent"
                and curr_gene not in genes_in_region
            ):
                right_border.append(curr_gene)
                single_copy_persistent_count += 1
            elif (
                curr_gene is not None
                and curr_gene not in genes_in_region
                and not return_only_persistents
            ):
                right_border.append(curr_gene)
            pos += 1
            if pos == self.contig.number_of_genes and self.contig.is_circular:
                pos = -1
            if pos == init:
                break  # looped around the contig

        border = [left_border, right_border]
        return border


class Spot(MetaFeatures):
    """
    The 'Spot' class represents a region of genomic plasticity.

    Methods:
        - 'regions': the property that generates the regions in the spot.
        - 'families': the property that generates the gene families in the spot.
        - 'spot_2_families': add to Gene Families a link to spot.
        - 'borders': Extracts all the borders of all RGPs belonging to the spot
        - 'get_uniq_to_rgp': Get dictionary with a representing RGP as key, and all identical RGPs as value
        - 'get_uniq_ordered_set': Get an Iterable of all the unique syntenies in the spot
        - 'get_uniq_content': Get an Iterable of all the unique rgp (in terms of gene family content) in the spot
        - 'count_uniq_content':  Get a counter of uniq RGP and number of identical RGP (in terms of gene family content)
        - 'count_uniq_ordered_set': Get a counter of uniq RGP and number of identical RGP (in terms of synteny content)

    Fields:
        - 'ID': Identifier of the spot
    """

    def __init__(self, spot_id: int):
        """Constructor method

        :param spot_id: Identifier of the spot
        """
        if not isinstance(spot_id, int):
            raise TypeError(
                f"Spot identifier must be an integer. Given type is {type(spot_id)}"
            )
        super().__init__()
        self.ID = spot_id
        self._region_getter = {}
        self._uniqOrderedSet = {}
        self._uniqContent = {}

    def __repr__(self) -> str:
        """Spot representation"""
        return f"Spot {self.ID} - #RGP: {len(self)}"

    def __str__(self):
        """String representation of the spot"""
        return f"spot_{self.ID}"

    def __setitem__(self, name: str, region: Region):
        """Set the region belonging to the spot

        :param name: Name of the region
        :param region: Region to add in the spot

        :raises KeyError: Name of the region is already in the spot for a different region
        """
        if name in self._region_getter and self[name] != region:
            raise KeyError("A Region with the same name already exist in spot")

        if not region.projected and region.spot is not None and region.spot != self:
            # In normal cases, a region should only belong to one spot. However, an exception arises in the projection command,
            # where a projected RGP might link two spots in the spot graph.
            # To handle this scenario without triggering failure, we check the 'projected' attribute of the given region.

            raise ValueError(
                f"The region '{region.name}' is already associated with spot '{region.spot.ID}' while being associated with spot '{self.ID}'. "
                "A region should only belong to one spot."
            )

        self._region_getter[name] = region
        region.spot = self

    def __getitem__(self, name) -> Region:
        """Get the region with the given name

        :param name: Name of the wanted region

        :return: Region in the spot for the given name

        :raises KeyError: Name does not exist in the spot
        :raises TypeError: Name is not a string
        """
        if not isinstance(name, str):
            raise TypeError(
                f"Name of the region must be a string. The provided type was {type(name)}"
            )
        try:
            region = self._region_getter[name]
        except KeyError:
            raise KeyError(f"Region with {name} does not exist in spot")
        else:
            return region

    def __delitem__(self, name):
        """Delete the region for the given name

        :param name: Name of the wanted region

        :raises KeyError: Name does not exist in the spot
        :raises TypeError: Name is not a string
        """
        if not isinstance(name, str):
            raise TypeError(
                f"Name of the region must be a string. The provided type was {type(name)}"
            )
        try:
            del self._region_getter[name]
        except KeyError:
            raise KeyError(f"Region with {name} does not exist in spot")

    def __len__(self) -> int:
        """Get the number of regions in the spot"""
        return len(self._region_getter)

    def add(self, region: Region):
        """Add a region to the spot.
        Alias more readable for setitem

        :param region: Region to add in the spot

        :raises TypeError: Region is not an instance Region
        """
        if not isinstance(region, Region):
            raise TypeError(
                f"A Region object is expected to be added to the spot. find type is {type(region)}"
            )
        self[region.name] = region

    def get(self, name: str) -> Region:
        """Get a region by its name.
        Alias more readable for getitem

        :param name: Name of the region

        :return: Wanted region
        """
        return self[name]

    def remove(self, name: str):
        """Remove a region by its name.
        Alias more readable for delitem

        :param name: Name of the region
        """
        del self[name]

    @property
    def regions(self) -> Generator[Region, None, None]:
        """Generates the regions in the spot

        :return: Regions in the spot
        """
        yield from self._region_getter.values()

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """Get the gene families in the RGP

        :return: Family in the spot
        """
        families = set()
        for region in self.regions:
            for family in region.families:
                if family not in families:
                    families.add(family)
                    yield family

    @property
    def number_of_families(self) -> int:
        """Get the number of different families in the spot

        :return: Number of families
        """
        return len({family for region in self.regions for family in region.families})

    def spot_2_families(self):
        """Add to Gene Families a link to spot"""
        for family in self.families:
            family.add_spot(self)

    def borders(
        self, set_size: int, multigenics
    ) -> List[List[int, List[GeneFamily], List[GeneFamily]]]:
        """Extracts all the borders of all RGPs belonging to the spot

        :param set_size: Number of genes to get
        :param multigenics: pangenome graph multigenic persistent families

        :return: Families that bordering spot
        """
        all_borders = [
            rgp.get_bordering_genes(set_size, multigenics) for rgp in self.regions
        ]

        family_borders = []
        for borders in all_borders:
            new = True
            curr_set = [
                [gene.family for gene in borders[0]],
                [gene.family for gene in borders[1]],
            ]
            for i, (c, former_borders) in enumerate(family_borders):
                if former_borders == curr_set or former_borders == curr_set[::-1]:
                    family_borders[i][0] += 1
                    new = False
                    break
            if new:
                family_borders.append([1, curr_set])

        return family_borders

    def _mk_uniq_ordered_set_obj(self):
        """cluster RGP into groups that have an identical synteny"""
        for rgp in self.regions:
            z = True
            for seen_rgp in self._uniqOrderedSet:
                if rgp == seen_rgp:
                    z = False
                    self._uniqOrderedSet[seen_rgp].add(rgp)
            if z:
                self._uniqOrderedSet[rgp] = {rgp}

    def _get_ordered_set(self) -> Dict[Region, Set[Region]]:
        """Creates the _uniqSyn object if it was never computed. Return it in any case

        :return: RGP groups that have an identical synteny
        """
        if len(self._uniqOrderedSet) == 0:
            self._mk_uniq_ordered_set_obj()
        return self._uniqOrderedSet

    def get_uniq_to_rgp(self) -> Dict[Region, Set[Region]]:
        """Get dictionary with a representing RGP as the key, and all identical RGPs as value

        :return: Dictionary with a representing RGP as the key, and set of identical RGPs as value
        """
        return self._get_ordered_set()

    def get_uniq_ordered_set(self) -> Set[Region]:
        """Get an Iterable of all the unique syntenies in the spot

        :return: Iterable of all the unique syntenies in the spot
        """
        return set(self._get_ordered_set().keys())

    def _mk_uniq_content(self):
        """cluster RGP into groups that have identical gene content"""
        for rgp in self.regions:
            z = True
            for seen_rgp in self._uniqContent:
                if rgp.families == seen_rgp.families:
                    z = False
                    self._uniqContent[seen_rgp].add(rgp)
            if z:
                self._uniqContent[rgp] = {rgp}

    def _get_content(self) -> Dict[Region, Set[Region]]:
        """Creates the _uniqContent object if it was never computed.

        :return: RGP groups that have identical gene content
        """
        if len(self._uniqContent) == 0:
            self._mk_uniq_content()
        return self._uniqContent

    def get_uniq_content(self) -> Set[Region]:
        """Get an Iterable of all the unique rgp (in terms of gene family content) in the spot

        :return: Iterable of all the unique rgp (in terms of gene family content) in the spot
        """
        return set(self._get_content().keys())

    def count_uniq_content(self) -> dict:
        """
        Get a counter of uniq RGP and number of identical RGP (in terms of gene family content)

        :return: Dictionary with a representative rgp as the key and number of identical rgp as value
        """
        return {key: len(val) for key, val in self._get_content().items()}

    def count_uniq_ordered_set(self):
        """
        Get a counter of uniq RGP and number of identical RGP (in terms of synteny content)

        :return: Dictionary with a representative rgp as the key and number of identical rgp as value
        """
        return {key: len(val) for key, val in self._get_ordered_set().items()}


class Module(MetaFeatures):
    """The `Module` class represents a module in a pangenome analysis.

    The `Module` class has the following attributes:
    - `ID`: An integer identifier for the module.
    - `bitarray`: A bitarray representing the presence/absence of the gene families in an organism.

    The `Module` class has the following methods:
    - `families`: Returns a generator that yields the gene families in the module.
    - `mk_bitarray`: Generates a bitarray representing the presence/absence of the gene families in an organism using the provided index.
    """

    def __init__(self, module_id: int, families: set = None):
        """Constructor method

        :param module_id: Module identifier
        :param families: Set of families which define the module
        """
        if not isinstance(module_id, int):
            raise TypeError(
                f"Module identifier must be an integer. Given type is {type(module_id)}"
            )
        super().__init__()
        self.ID = module_id
        self._families_getter = {}
        self.bitarray = None
        if families is not None:
            for family in families:
                self.add(family)

    def __repr__(self) -> str:
        """Module representation"""
        return f"Module {self.ID} - #Families: {len(self)}"

    def __str__(self) -> str:
        """String representation of the module"""
        return f"module_{self.ID}"

    def __hash__(self) -> int:
        """Create a hash value for the module"""
        return id(self)

    def __len__(self) -> int:
        """Get the number of families in the module"""
        return len(self._families_getter)

    def __eq__(self, other: Module) -> bool:
        """
        Test whether two Module objects have the same gene families

        :param other: Another module to test equality

        :return: Equal or not

        :raises TypeError: Try to compare a module with another type object
        """
        if not isinstance(other, Module):
            raise TypeError(
                f"Another module is expected to be compared to the first one. You give a {type(other)}"
            )
        return set(self.families) == set(other.families)

    def __setitem__(self, name: str, family: GeneFamily):
        """Set a gene family in the module

        :param name: Name of the family
        :param family: Gene family belonging to the module

        :raises TypeError: Family is not instance GeneFamily
        :raises KeyError: Another family with the same name already exists in the module
        """
        if name in self._families_getter and self[name] != family:
            raise KeyError(
                "A different gene family with the same name already exist in the module"
            )
        self._families_getter[name] = family
        family.set_module(self)

    def __getitem__(self, name) -> GeneFamily:
        """Get the gene family for the given name in the module

        :param name: Name of the gene family

        :return: Gene family with the given name

        :raises KeyError: Family with the given name does not exist in the module
        """
        try:
            family = self._families_getter[name]
        except KeyError:
            raise KeyError(
                f"There isn't gene family with the name {name} in the module"
            )
        else:
            return family

    def __delitem__(self, name):
        """Remove the gene family for the given name in the module

        :param name: Name of the gene family

        :raises KeyError: Family with the given name does not exist in the module
        """
        try:
            fam = self._families_getter[name]
        except KeyError:
            raise KeyError(
                f"There isn't gene family with the name {name} in the module"
            )
        else:
            del self._families_getter[name]
            fam._module = None  # TODO define method to remove a module from family

    def add(self, family: GeneFamily):
        """Add a family to the module.
        Alias more readable for setitem

        :param family: Region to add in the spot

        :raises TypeError: Region is not an instance Region
        """
        if not isinstance(family, GeneFamily):
            raise TypeError(
                f"A gene family is expected to be added to module. Given type was {type(family)}"
            )
        self[family.name] = family

    def get(self, name: str) -> GeneFamily:
        """Get a family by its name.
        Alias more readable for getitem

        :param name: Name of the family

        :return: Wanted family
        """
        return self[name]

    def remove(self, name: str):
        """Remove a family by its name.
        Alias more readable for delitem

        :param name: Name of the family
        """
        del self[name]

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """Generator of the family in the module

        :return: Families belonging to the module
        """
        yield from self._families_getter.values()

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """Returns all the Organisms that have this module

        :return: Organisms that have this module
        """
        organisms = set()
        for fam in self.families:
            organisms |= set(fam.organisms)
        yield from organisms

    def mk_bitarray(self, index: Dict[GeneFamily, int], partition: str = "all"):
        """Produces a bitarray representing the presence / absence of families in the organism using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.

        :param partition: filter module by partition
        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`
        """
        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == "all":
            logging.getLogger("PPanGGOLiN").debug("all")
            for fam in self.families:
                self.bitarray[index[fam]] = 1
        elif partition == "persistent":
            logging.getLogger("PPanGGOLiN").debug("persistent")
            for fam in self.families:
                if fam.named_partition in ["persistent"]:
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
            raise Exception(
                "There is not any partition corresponding please report a github issue"
            )


class GeneContext:
    """
    Represent a gene context which is a collection of gene families related to a specific genomic context.

    Methods
    - families: Generator that yields all the gene families in the gene context.
    - add_context_graph: Add a context graph corresponding to the gene context.
    - add_family: Add a gene family to the gene context.

    Fields
    - gc_id: The identifier of the gene context.
    - graph: context graph corresponding to the gene context
    """

    def __init__(
        self,
        gc_id: int,
        families: Set[GeneFamily] = None,
        families_of_interest: Set[GeneFamily] = None,
    ):
        """Constructor method

        :param gc_id: Identifier of the gene context.
        :param families: Gene families related to the gene context.
        :param families_of_interest: Input families for which the context is being searched.
        """

        if not isinstance(gc_id, int):
            raise TypeError(
                f"Gene context identifier must be an integer. Given type is {type(gc_id)}"
            )

        self.ID = gc_id
        self._families_getter = {}
        self.families_of_interest = families_of_interest
        self._graph = None
        if families is not None:
            if not all(isinstance(fam, GeneFamily) for fam in families):
                raise Exception(
                    "You provided elements that were not GeneFamily objects. "
                    "GeneContexts are only made of GeneFamily objects."
                )
            self._families_getter = {family.name: family for family in families}

    def __repr__(self) -> str:
        """Context representation"""
        return f"Context {self.ID} - #Families: {len(self)}"

    def __str__(self) -> str:
        """String representation of the gene context"""
        return f"GC_{str(self.ID)}"

    def __hash__(self) -> int:
        """Create a hash value for the region"""
        return id(self)

    def __len__(self) -> int:
        """Get the number of families in the context"""
        return len(self._families_getter)

    def __eq__(self, other: GeneContext) -> bool:
        """
        Test whether two gene context objects have the same gene families

        :param other: Another gene context to test equality

        :return: Equal or not

        :raises TypeError: Try to compare a gene context with another type object
        """
        if not isinstance(other, GeneContext):
            raise TypeError(
                f"Another context is expected to be compared to the first one. You give a {type(other)}"
            )
        return set(self.families) == set(other.families)

    def __setitem__(self, name, family):
        """Set a gene family in the gene context

        :param name: Name of the family
        :param family: Gene family belonging to the context

        :raises TypeError: Family is not instance GeneFamily
        :raises KeyError: Another family with the same name already exists in the context
        """
        if not isinstance(family, GeneFamily):
            raise TypeError(
                f"A gene family is expected to be added to gene context. Given type was {type(family)}"
            )
        if name in self._families_getter and self[name] != family:
            raise KeyError(
                "A different gene family with the same name already exist in the gene context"
            )
        self._families_getter[name] = family

    def __getitem__(self, name) -> GeneFamily:
        """Get the gene family for the given name in the context

        :param name: Name of the gene family

        :return: Gene family with the given name

        :raises KeyError: Family with the given name does not exist in the context
        """
        try:
            family = self._families_getter[name]
        except KeyError:
            raise KeyError(
                f"There isn't gene family with the name {name} in the gene context"
            )
        else:
            return family

    def __delitem__(self, name):
        """Remove the gene family for the given name in the context

        :param name: Name of the gene family

        :raises KeyError: Family with the given name does not exist in the context
        """
        try:
            del self._families_getter[name]
        except KeyError:
            raise KeyError(
                f"There isn't gene family with the name {name} in the gene context"
            )

    @property
    def graph(self):
        if self._graph is None:
            raise ValueError("Graph has not been added to the context")
        return self._graph

    @graph.setter
    def graph(self, graph: nx.Graph):
        """
        Add a context graph to the gene context.

        :param graph: The context graph.
        """
        if not isinstance(graph, nx.Graph):
            logging.getLogger("PPanGGOLiN").debug(f"given type: {type(graph)}")
            raise TypeError("Context graph must be a networkx graph object.")
        self._graph = graph

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """Generator of the family in the context

        :return: Gene families belonging to the context
        """
        yield from self._families_getter.values()

    def add_family(self, family: GeneFamily):
        """
        Add a gene family to the gene context.

        :param family: The gene family to add.
        """
        if not isinstance(family, GeneFamily):
            raise Exception(
                "You did not provide a GeneFamily object. "
                "GeneContexts are only made of GeneFamily objects."
            )
        self[family.name] = family
