#!/usr/bin/env python3
# coding: utf8

# default libraries
from __future__ import annotations
import logging
from collections.abc import Iterable

# installed libraries
from typing import Dict, Set

import gmpy2

# local libraries
from ppanggolin.genome import Gene, Organism, Contig
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.metadata import MetaFeatures


class Region(MetaFeatures):
    """
    This class represent a region of genomic plasticity.

    :param region_id: identifier of the region
    """

    def __init__(self, region_id: str):
        super().__init__()
        self.genes = []
        self.name = region_id
        self.score = 0

    def __str__(self):
        return self.name

    def __hash__(self):
        return id(self)

    def __eq__(self, other: Region) -> bool:
        """
        Expects another Region type object. Will test whether two Region objects have the same gene families

        :param other: Other region to test equality of region

        :return: equal or not
        """
        if not isinstance(other, Region):
            raise TypeError(f"'Region' type object was expected, but '{type(other)}' type object was provided.")
        if [gene.family for gene in self.genes] == [gene.family for gene in other.genes]:
            return True
        if [gene.family for gene in self.genes] == [gene.family for gene in other.genes[::-1]]:
            return True
        return False

    def __len__(self):
        return len(self.genes)

    def __getitem__(self, index):
        return self.genes[index]

    def append(self, gene: Gene):
        """allowing only gene-class objects in a region

        :param gene: gene which will be added

        :raise TypeError: If gene is not Gene type raise TypeError
        """

        if isinstance(gene, Gene):
            self.genes.append(gene)
            gene.RGP.add(self)
        else:
            raise TypeError(f"Unexpected class / type for {type(gene)} when adding it to a RGP")

    @property
    def families(self) -> Set[GeneFamily]:
        """Get the gene families in the RGP

        :return: Set of gene families
        """
        return {gene.family for gene in self.genes}

    @property
    def start(self) -> int:
        """ Get RGP starting position

        :return: Start position
        """
        return min(self.genes, key=lambda x: x.start).start

    @property  # TODO try to change start with this method
    def start_gene(self) -> Gene:
        """ Get RGP starting gene

        :return: Start gene
        """
        return min(self.genes, key=lambda x: x.position)

    @property
    def stop_gene(self) -> Gene:
        """ Get RGP stoping position

        :return: Stoping position
        """
        return max(self.genes, key=lambda x: x.position)

    @property
    def stop(self):
        """ Get RGP stoping position

        :return: Stop position
        """
        return max(self.genes, key=lambda x: x.stop).stop

    @property
    def organism(self) -> Organism:
        """ Get the Organism link to RGP

        :return: Organism
        """
        return self.genes[0].organism

    @property
    def contig(self) -> Contig:
        """ Get the Contig link to RGP

        :return: Contig
        """
        return self.genes[0].contig

    @property
    def is_whole_contig(self) -> bool:
        """Indicates if the region is an entire contig

        :return: True if whole contig
        """
        if self.start_gene.position == 0 and self.stop_gene.position == len(self.contig.genes) - 1:
            return True
        return False

    @property
    def is_contig_border(self) -> bool:
        """Indicates if the region is bordering a contig

        :return: True if bordering
        """
        if len(self.genes) == 0:
            raise Exception("Your region has no genes. Something wrong happenned.")
        if (self.start_gene.position == 0 and not self.contig.is_circular) or \
                (self.stop_gene.position == len(self.contig.genes) - 1 and not self.contig.is_circular):
            return True
        return False

    def get_rnas(self) -> set:
        """ Get RNA in region

        :return: Set of RNA
        """
        rnas = set()
        for rna in self.contig.RNAs:
            if self.start < rna.start < self.stop:
                rnas.add(rna)
        return rnas

    def get_bordering_genes(self, n: int, multigenics: set) -> list:
        """ Get the bordered genes in the region

        :param n: number of genes to get
        :param multigenics: pangenome graph multigenic persistent families

        :return: A list of bordering gene in start and stop position List[List[Start Gene], [Stop Gene]]
        """
        border = [[], []]
        pos = self.start_gene.position
        init = pos
        while len(border[0]) < n and (pos != 0 or self.contig.is_circular):
            curr_gene = None
            if pos == 0:
                if self.contig.is_circular:
                    curr_gene = self.contig.genes[-1]
            else:
                curr_gene = self.contig.genes[pos - 1]
            if curr_gene is not None and curr_gene.family not in multigenics and \
                    curr_gene.family.named_partition == "persistent":
                border[0].append(curr_gene)
            pos -= 1
            if pos == -1 and self.contig.is_circular:
                pos = len(self.contig.genes)
            if pos == init:
                break  # looped around the contig
        pos = self.stop_gene.position
        init = pos
        while len(border[1]) < n and (pos != len(self.contig.genes) - 1 or self.contig.is_circular):
            curr_gene = None
            if pos == len(self.contig.genes) - 1:
                if self.contig.is_circular:
                    curr_gene = self.contig.genes[0]
            else:
                curr_gene = self.contig.genes[pos + 1]
            if curr_gene is not None and curr_gene.family not in multigenics:
                border[1].append(curr_gene)
            pos += 1
            if pos == len(self.contig.genes) and self.contig.is_circular:
                pos = -1
            if pos == init:
                break  # looped around the contig
        return border


class Spot(MetaFeatures):
    """
    This class represent a hotspot.

    :param spot_id: identifier of the spot
    """
    def __init__(self, spot_id):
        super().__init__()
        self.ID = spot_id
        self.regions = set()
        self._uniqOrderedSet = {}
        self._compOrderedSet = False
        self._uniqContent = {}
        self._compContent = False

    def __str__(self):
        return f'spot_{str(self.ID)}'

    @property
    def families(self) -> set:
        """Get the gene families in the RGP

        :return: Set of gene families
        """

        union = set()
        for region in self.regions:
            union |= region.families
        return union

    def add_regions(self, regions):
        """
        Adds region(s) contained in an Iterable to the spot which all have the same bordering persistent genes
        provided with 'borders'

        :param regions: Iterable list of RGP to add to spot
        """
        if isinstance(regions, Iterable):
            for region in regions:
                self.add_region(region)
        else:
            raise Exception("The provided 'regions' variable was not an Iterable")

    def add_region(self, region):
        """
        Add one RGP to the spot

        :param region: RGP to add to spot
        """
        if isinstance(region, Region):
            self.regions.add(region)

    def spot_2_families(self):
        """Add to Gene Families a link to spot"""
        for family in self.families:
            family.spot.add(self)

    def borders(self, set_size: int, multigenics):
        """ Extracts all the borders of all RGPs belonging to the spot

        :param set_size: number of genes to get
        :param multigenics: pangenome graph multigenic persistent families

        :return: families that bordering spot
        """
        all_borders = []
        for rgp in self.regions:
            all_borders.append(rgp.get_bordering_genes(set_size, multigenics))

        family_borders = []
        for borders in all_borders:
            new = True
            curr_set = [[gene.family for gene in borders[0]], [gene.family for gene in borders[1]]]
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

    def _get_content(self):
        """Creates the _uniqContent object if it was never computed. Return it in any case

        :return: RGP groups that have identical gene content
        """
        if not self._compContent:
            self._mk_uniq_content()
            self._compContent = True
        return self._uniqContent

    def _get_ordered_set(self):
        """ Creates the _uniqSyn object if it was never computed. Return it in any case

        :return: RGP groups that have an identical synteny
        """
        if not self._compOrderedSet:
            self._mk_uniq_ordered_set_obj()
            self._compOrderedSet = True
        return self._uniqOrderedSet

    def get_uniq_to_rgp(self) -> dict:
        """ Get dictionnary with a representing RGP as key, and all identical RGPs as value

        :return: Dictionnary with a representing RGP as key, and all identical RGPs as value
        """
        return self._get_ordered_set()

    def get_uniq_ordered_set(self):
        """Get an Iterable of all the unique syntenies in the spot

        :return: Iterable of all the unique syntenies in the spot
        """
        return set(self._get_ordered_set().keys())

    def get_uniq_content(self):
        """ Get an Iterable of all the unique rgp (in terms of gene family content) in the spot

        :return: Iterable of all the unique rgp (in terms of gene family content) in the spot
        """
        return set(self._get_content().keys())

    def count_uniq_content(self) -> dict:
        """
        Get a counter of uniq RGP and number of identical RGP (in terms of gene family content)

        :return: dictionary with a representative rgp as key and number of identical rgp as value
        """
        return dict([(key, len(val)) for key, val in self._get_content().items()])

    def count_uniq_ordered_set(self):
        """
        Get a counter of uniq RGP and number of identical RGP (in terms of synteny content)

        :return: dictionary with a representative rgp as key and number of identical rgp as value
        """
        return dict([(key, len(val)) for key, val in self._get_ordered_set().items()])


class Module(MetaFeatures):
    """
    This class represent a hotspot.

    :param module_id: identifier of the module
    :param families: Set of families which define the module
    """
    def __init__(self, module_id: int, families: set = None):
        """
        'core' are gene families that define the module.
        'associated_families' are gene families that you believe are associated to the module in some way,
        but do not define it.
        """
        super().__init__()
        self.ID = module_id
        self._families = set()
        if families is not None:
            if not all(isinstance(fam, GeneFamily) for fam in families):
                raise Exception("You provided elements that were not GeneFamily object. "
                                "Modules are only made of GeneFamily")
            self._families |= set(families)
        self.bitarray = None

    @property
    def families(self) -> Set[GeneFamily]:
        # TODO made as generator
        return self._families

    def __str__(self):
        return f'module_{str(self.ID)}'

    def add_family(self, family: GeneFamily):
        """
        Add a family to the module

        :param family: the family that will ba added to the module
        """
        if not isinstance(family, GeneFamily):
            raise Exception("You did not provide a GenFamily object. Modules are only made of GeneFamily")
        family.modules.add(self)
        self._families.add(family)

    def mk_bitarray(self, index: Dict[Organism, int], partition: str = 'all'):
        """Produces a bitarray representing the presence / absence of families in the organism using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.

        :param partition: filter module by partition
        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`
        """
        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == 'all':
            logging.getLogger("PPanGGOLiN").debug("all")
            for fam in self.families:
                self.bitarray[index[fam]] = 1
        elif partition == 'persistent':
            logging.getLogger("PPanGGOLiN").debug("persistent")
            for fam in self.families:
                if fam.named_partition in ['persistent']:
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


class GeneContext:
    """
    A class used to represent a gene context

    :param gc_id : identifier of the Gene context
    :param families: Gene families related to the GeneContext
    """

    def __init__(self, gc_id: int, families: set = None):
        self.ID = gc_id
        self.families = set()
        if families is not None:
            if not all(isinstance(fam, GeneFamily) for fam in families):
                raise Exception("You provided elements that were not GeneFamily object."
                                " GeneContext are only made of GeneFamily")
            self.families |= set(families)

    def __str__(self):
        return f'GC_{str(self.ID)}'

    def add_family(self, family: GeneFamily):
        """
        Allow to add one family in the GeneContext
        :param family: family to add
        """
        if not isinstance(family, GeneFamily):
            raise Exception("You did not provide a GenFamily object. Modules are only made of GeneFamily")
        self.families.add(family)
