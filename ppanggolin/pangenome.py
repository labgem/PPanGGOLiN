#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Iterator, List, Union, Dict, Set, Iterable
from pathlib import Path

# local libraries
from ppanggolin.genome import Organism, Gene
from ppanggolin.region import Region, Spot, Module
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.edge import Edge


class Pangenome:
    """
    This is a class representing your pangenome. It is used as a basic unit for all the analysis to access to the
    different elements of your pangenome, such as organisms, contigs, genes or gene families. It has setter and getter
    methods for most elements in your pangenome, and you can use those to add new elements to it,
    or get objects that have a specific identifier to manipulate them directly.
    """

    def __init__(self):
        """Constructor method.
        """
        self.file = None

        # basic parameters
        self._famGetter = {}
        self._org_index = None
        self._fam_index = None
        self.max_fam_id = 0
        self._orgGetter = {}
        self._edgeGetter = {}
        self._regionGetter = {}
        self._spots = set()
        self._modules = set()

        self.status = {
            'genomesAnnotated': "No",
            'geneSequences': "No",
            'genesClustered': "No",
            'defragmented': "No",
            'geneFamilySequences': "No",
            'neighborsGraph': "No",
            'partitioned': "No",
            'predictedRGP': "No",
            'spots': "No",
            'modules': 'No'
        }
        self.parameters = {}

    def add_file(self, pangenome_file: Path):
        """Links an HDF5 file to the pangenome. If needed elements will be loaded from this file,
        and anything that is computed will be saved to this file when
        :func:`ppanggolin.formats.writeBinaries.writePangenome` is called.

        :param pangenome_file: A string representing filepath to hdf5 pangenome file to be either used or created
        """
        from ppanggolin.formats.readBinaries import get_status
        # importing on call instead of importing on top to avoid cross-reference problems.
        get_status(self, pangenome_file)
        self.file = pangenome_file.absolute().as_posix()

    """ Gene Methods"""

    @property
    def genes(self) -> list:
        """Creates the geneGetter if it does not exist, and returns all the genes of all organisms in the pangenome.
        
        :return: list of :class:`ppanggolin.genome.Gene`
        """
        try:
            return list(self._geneGetter.values())
        except AttributeError:  # in that case the gene getter has not been computed
            self._mk_gene_getter()  # make it
            return self.genes  # return what was expected

    def _yield_genes(self) -> Iterator[Gene]:
        """ Use a generator to get all the genes of a pangenome

        :return: an iterator of Gene
        """
        if self.number_of_organisms() > 0:  # if we have organisms, they're supposed to have genes
            for org in self.organisms:
                for contig in org.contigs:
                    for gene in contig.genes:
                        yield gene
        elif self.number_of_gene_families() > 0:
            # we might have no organism loaded, in that case there are gene families.
            for geneFam in self.gene_families:
                for gene in geneFam.genes:
                    yield gene

    def _mk_gene_getter(self):
        """
        Builds the attribute _geneGetter of the pangenome

        Since the genes are never explicitly 'added' to a pangenome (but rather to a gene family, or a contig),
        the pangenome cannot directly extract a gene from a geneID since it does not 'know' them.
        if at some point we want to extract genes from a pangenome we'll create a geneGetter.
        The assumption behind this is that the pangenome has been filled and no more gene will be added.
        """
        self._geneGetter = {}
        for gene in self._yield_genes():
            self._geneGetter[gene.ID] = gene

    def get_gene(self, gene_id: str) -> Gene:
        """returns the gene that has the given geneID

        :param gene_id: The gene ID to look for

        :return: returns the gene that has the ID `geneID`

        :raises KeyError: If the `geneID` is not in the pangenome
        """
        try:
            return self._geneGetter[gene_id]
        except AttributeError:
            # in that case, either the gene getter has not been computed, or the geneID is not in the pangenome.
            self._mk_gene_getter()  # make it
            return self.get_gene(gene_id)  # return what was expected. If geneID does not exist it will raise an error.
        except KeyError:
            raise KeyError(f"{gene_id} does not exist in the pangenome.")

    def number_of_gene(self) -> int:
        """Returns the number of gene present in the pangenome

        :return: the number of gene families
        """
        try:
            return len(self._geneGetter)
        except AttributeError:  # in that case the gene getter has not been computed
            self._mk_gene_getter()  # make it
            return len(self._geneGetter)

    """Gene families methods"""

    @property
    def gene_families(self) -> List[GeneFamily]:
        """returns all the gene families in the pangenome
        
        :return: list of :class:`ppanggolin.geneFamily.GeneFamily`
        """
        return list(self._famGetter.values())

    def _create_gene_family(self, name: str) -> GeneFamily:
        """Creates a gene family object with the given `name`

        :param name: the name to give to the gene family. Must not exist already.

        :return: the created GeneFamily object
        """
        new_fam = GeneFamily(family_id=self.max_fam_id, name=name)
        self.max_fam_id += 1
        self._famGetter[new_fam.name] = new_fam
        return new_fam

    def number_of_gene_families(self) -> int:
        """Returns the number of gene families present in the pangenome

        :return: the number of gene families
        """
        return len(self._famGetter)

    def get_gene_family(self, name: str) -> GeneFamily:
        """returns the gene family that has the given `name`

        :param name: The gene family name to look for

        :return: returns the gene family that has the name `name`
        """
        return self._famGetter[name]

    def add_gene_family(self, name: str):
        """
        Get the :class:`ppanggolin.geneFamily.GeneFamily` object that has the given `name`.
        If it does not exist, creates it.

        :param name: The gene family name to get if it exists, and create otherwise.

        :return: GeneFamily object.
        """
        fam = self._famGetter.get(name)
        if fam is None:
            fam = self._create_gene_family(name)
        return fam

    """Graph methods"""

    @property
    def edges(self) -> list:
        """returns all the edges in the pangenome graph
        
        :return: list of :class:`ppanggolin.pangenome.Edge`
        """
        return list(self._edgeGetter.values())

    def add_edge(self, gene1: Gene, gene2: Gene) -> Edge:
        """
        Adds an edge between the two gene families that the two given genes belong to.
        Genes object are expected, and they are also expected to have a family assigned

        :param gene1: The first gene
        :param gene2: The second gene

        :return: the created Edge
        """
        key = frozenset([gene1.family, gene2.family])
        edge = self._edgeGetter.get(key)
        if edge is None:
            edge = Edge(gene1, gene2)
            self._edgeGetter[key] = edge
        else:
            edge.add_genes(gene1, gene2)
        return edge

    def number_of_edge(self) -> int:
        """Returns the number of edge present in the pangenome

        :return: the number of gene families
        """
        return len(self._edgeGetter)

    """Organism methods"""

    @property
    def organisms(self) -> List[Organism]:
        """returns all the organisms in the pangenome
        
        :return: list of :class:`ppanggolin.genome.Organism`
        """
        return list(self._orgGetter.values())

    def number_of_organisms(self) -> int:
        """Returns the number of organisms present in the pangenome
        
        :return: the number of organism
        """
        return len(self._orgGetter)

    def get_organism(self, org_name: str) -> Organism:
        """
        Get an organism that is expected to be in the pangenome using its name, which is supposedly unique.
        Raises an error if the organism does not exist.

        :param org_name: Name of the Organism to get

        :return: The related Organism object

        :raises KeyError: If the provided name is not in the pangenome
        """
        try:
            return self._orgGetter[org_name]
        except KeyError:
            raise KeyError(f"{org_name} does not seem to be in your pangenome")

    def add_organism(self, new_org: Union[Organism, str]) -> Organism:
        """
        adds an organism that did not exist previously in the pangenome if an Organism object is provided.
        If an organism with the same name exists it will raise an error.
        If a str object is provided, will return the corresponding organism that has this name
        OR create a new one if it does not exist.

        :param new_org: Organism to add to the pangenome

        :return: The created organism

        :raises TypeError: if the provided `newOrg` is neither a str nor a :class:`ppanggolin.genome.Organism`
        """
        if isinstance(new_org, Organism):
            old_len = len(self._orgGetter)
            self._orgGetter[new_org.name] = new_org
            if len(self._orgGetter) == old_len:
                raise KeyError(f"Redondant organism name was found ({new_org.name})."
                               f"All of your organisms must have unique names.")
        elif isinstance(new_org, str):
            org = self._orgGetter.get(new_org)
            if org is None:
                org = Organism(new_org)
                self._orgGetter[org.name] = org
            new_org = org
        else:
            raise TypeError("Provide an Organism object or a str that will serve as organism name")
        return new_org

    def get_org_index(self) -> Dict[Organism, int]:  # will not make a new index if it exists already
        """Creates an index for Organisms (each organism is assigned an Integer).

        :return: The index of organisms in pangenome
        """
        if self._org_index is None:  # then the index does not exist yet
            self._org_index = {}
            for index, org in enumerate(self.organisms):
                self._org_index[org] = index
        return self._org_index

    def compute_family_bitarrays(self, part: str = 'all') -> Dict[Organism, int]:
        """
        Based on the index generated by get_org_index, generate a bitarray for each gene family.
        If the family j is present in the organism with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :param part: Filter the organism in function of the given partition

        :return: the index of organisms in pangenome
        """
        if self._org_index is None:
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.get_org_index()
        for fam in self.gene_families:
            fam.mk_bitarray(self._org_index, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._org_index

    def get_fam_index(self) -> Dict[GeneFamily, int]:  # will not make a new index if it exists already
        """Creates an index for gene families (each family is assigned an Integer).

        :return: The index of families in pangenome
        """
        if self._fam_index is None:  # then the index does not exist yet
            self._fam_index = {}
            for index, fam in enumerate(self.gene_families):
                self._fam_index[fam] = index
        return self._fam_index

    def compute_org_bitarrays(self, part='all') -> Dict[GeneFamily, int]:
        """
        Based on the index generated by get_fam_index, generate a bitarray for each gene family.
        If the family j is present in the organism with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :param part: Filter the organism in function of the given partition

        :return: The index of gene families in pangenome
        """
        if self._fam_index is None:
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.get_fam_index()
        for org in self.organisms:
            org.mk_bitarray(index=self._fam_index, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._fam_index

    """RGP methods"""

    @property
    def regions(self) -> list:
        """returns all the regions (RGP) in the pangenome

        :return: list of RGP
        """
        return list(self._regionGetter.values())

    def get_region(self, region_name: str) -> Region:
        """Returns a region with the given region_name. Creates it if it does not exist.

        :param region_name: The name of the region to return

        :return: The region
        """
        try:
            return self._regionGetter[region_name]
        except KeyError:  # then the region is not stored in this pangenome.
            new_region = Region(region_name)
            self._regionGetter[region_name] = new_region
            return new_region

    def get_multigenics(self, dup_margin: float, persistent: bool = True) -> Set[GeneFamily]:
        """
        Returns the multigenic persistent families of the pangenome graph. A family will be considered multigenic
        if it is duplicated in more than `dup_margin` of the genomes where it is present.

        :param dup_margin: the ratio of presence in multicopy above which a gene family is considered multigenic
        :param persistent: if we consider only the persistent genes

        :return: set of gene families considered multigenic
        """
        multigenics = set()
        for fam in self.gene_families:
            if fam.named_partition == "persistent" or not persistent:
                dup = len([genes for org, genes in fam.get_org_dict().items() if
                           len([gene for gene in genes if not gene.is_fragment]) > 1])
                if (dup / len(fam.organisms)) >= dup_margin:  # tot / nborgs >= 1.05
                    multigenics.add(fam)
        # logging.info(f"{len(multigenics)} gene families are defined as being multigenic.
        # (duplicated in more than {dup_margin} of the genomes)")
        return multigenics

    def add_regions(self, region_group: Union[Region, Iterable[Region]]):
        """Takes an Iterable or a Region object and adds it to the pangenome

        :param region_group: a region or an Iterable of regions to add to the pangenome

        :raises TypeError: if regionGroup is neither a Region nor an Iterable[Region]
        """
        old_len = len(self._regionGetter)
        if isinstance(region_group, Iterable):
            for region in region_group:
                self._regionGetter[region.name] = region
            if len(self._regionGetter) != len(region_group) + old_len:
                raise Exception("Two regions had an identical name, which was unexpected.")
        elif isinstance(region_group, Region):
            self._regionGetter[region_group.name] = region_group
        else:
            raise TypeError(f"An iterable or a 'Region' type object were expected, "
                            f"but you provided a {type(region_group)} type object")

    def number_of_rgp(self) -> int:
        """Returns the number of gene families present in the pangenome

        :return: the number of gene families
        """
        return len(self._regionGetter)

    """Spot methods"""

    @property
    def spots(self) -> Set[Spot]:
        # TODO made as generator
        return self._spots

    def add_spots(self, spots: Iterable[Spot]):
        """Adds the given iterable of spots to the pangenome.

        :param spots: An iterable of :class:`ppanggolin.region.Spot`.
        """
        self._spots |= set(spots)

    def number_of_spots(self) -> int:
        """Returns the number of gene families present in the pangenome

        :return: the number of gene families
        """
        return len(self.spots)

    """Modules methods"""

    @property
    def modules(self) -> Set[Module]:
        # TODO made as generator
        return self._modules

    def add_modules(self, modules: Iterable[Module]):
        """Adds the given iterable of modules to the pangenome

        :param modules: an iterable of :class:`ppanggolin.module.Module`
        """
        self._modules |= set(modules)

    def compute_mod_bitarrays(self, part: str = 'all') -> Dict[GeneFamily, int]:
        """Based on the index generated by get_fam_index, generated a bitarray
        for each gene family present in modules.
        If the family j is present in the module with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :param part: Filter the organism in function of the given partition

        :return: A dictionary with Organism as key and int as value.
        """
        if self._fam_index is None:
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.get_fam_index()
        for module in self.modules:
            module.mk_bitarray(index=self._fam_index, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._fam_index

    def number_of_modules(self) -> int:
        """Returns the number of modules present in the pangenome

        :return: the number of modules
        """
        return len(self.modules)
