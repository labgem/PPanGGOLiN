#!/usr/bin/env python3
# coding: utf8

# default libraries
from collections.abc import Iterable

# local libraries
from ppanggolin.genome import Organism
from ppanggolin.region import Region
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.edge import Edge


class Pangenome:
    """
    This is a class representing your pangenome. It is used as a basic unit for all of the analysis to access to the
    different elements of your pangenome, such as organisms, contigs, genes or gene families. It has setter and getter
    methods for most elements in your pangenome and you can use those to add new elements to it,
    or get objects that have a specific identifier to manipulate them directly.
    """

    def __init__(self):
        """Constructor method.
        """
        self.file = None

        # basic parameters
        self._famGetter = {}
        self.max_fam_id = 0
        self._orgGetter = {}
        self._edgeGetter = {}
        self._regionGetter = {}
        self.spots = set()
        self.modules = set()

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

    def addFile(self, pangenomeFile):
        """Links an HDF5 file to the pangenome. If needed elements will be loaded from this file,
        and anything that is computed will be saved to this file when
        :func:`ppanggolin.formats.writeBinaries.writePangenome` is called.

        :param pangenomeFile: A string representing the filepath to the hdf5 pangenome file to be either used or created
        :type pangenomeFile: str
        """
        from ppanggolin.formats import \
            getStatus  # importing on call instead of importing on top to avoid cross-reference problems.
        getStatus(self, pangenomeFile)
        self.file = pangenomeFile

    """ Gene Methods"""
    @property
    def genes(self):
        """Creates the geneGetter if it does not exist, and returns all the genes of all organisms in the pangenome.
        
        :return: list of :class:`ppanggolin.genome.Gene`
        :rtype: list
        """
        try:
            return list(self._geneGetter.values())
        except AttributeError:  # in that case the gene getter has not been computed
            self._mkgeneGetter()  # make it
            return self.genes  # return what was expected

    def _yield_genes(self):
        """
        Use a generator to get all the genes of a pangenome

        :return: an iterator of :class:`ppanggolin.genome.Gene`
        :rtype: Iterator[:class:`ppanggolin.genome.Gene`]
        """
        if self.number_of_organisms() > 0:  # if we have organisms, they're supposed to have genes
            for org in self.organisms:
                for contig in org.contigs:
                    for gene in contig.genes:
                        yield gene
        elif self.number_of_geneFamilies() > 0:
            # we might have no organism loaded, in that case there are gene families.
            for geneFam in self.geneFamilies:
                for gene in geneFam.genes:
                    yield gene

    def _mkgeneGetter(self):
        """
            Builds the :attr:`ppanggolin.pangenome.Pangenome._geneGetter` of the pangenome

            Since the genes are never explicitly 'added' to a pangenome (but rather to a gene family, or a contig),
            the pangenome cannot directly extract a gene from a geneID since it does not 'know' them.
            if at some point we want to extract genes from a pangenome we'll create a geneGetter.
            The assumption behind this is that the pangenome has been filled and no more gene will be added.
        """
        self._geneGetter = {}
        for gene in self._yield_genes():
            self._geneGetter[gene.ID] = gene

    def getGene(self, geneID):
        """returns the gene that has the given `geneID`

        :param geneID: The gene ID to look for
        :type geneID: any
        :return: returns the gene that has the ID `geneID`
        :rtype: :class:`ppanggolin.genome.Gene`
        :raises KeyError: If the `geneID` is not in the pangenome
        """
        try:
            return self._geneGetter[geneID]
        except AttributeError:
            # in that case, either the gene getter has not been computed, or the geneID is not in the pangenome.
            self._mkgeneGetter()  # make it
            return self.getGene(
                geneID)  # return what was expected. If the geneID does not exist it will raise an error.
        except KeyError:
            raise KeyError(f"{geneID} does not exist in the pangenome.")

    """Gene families methods"""
    @property
    def geneFamilies(self):
        """returns all the gene families in the pangenome
        
        :return: list of :class:`ppanggolin.geneFamily.GeneFamily`
        :rtype: list
        """
        return list(self._famGetter.values())

    def _createGeneFamily(self, name):
        """Creates a gene family object with the given `name`

        :param name: the name to give to the gene family. Must not exist already.
        :type name: any
        :return: the created GeneFamily object
        :rtype: :class:`ppanggolin.geneFamily.GeneFamily`
        """
        newFam = GeneFamily(ID=self.max_fam_id, name=name)
        self.max_fam_id += 1
        self._famGetter[newFam.name] = newFam
        return newFam

    def number_of_geneFamilies(self):
        """Returns the number of gene families present in the pangenome

        :return: the number of gene families
        :rtype: int
        """
        return len(self._famGetter)

    def getGeneFamily(self, name):
        """returns the gene family that has the given `name`

        :param name: The gene family name to look for
        :type name: any
        :return: returns the gene family that has the name `name`
        :rtype: :class:`ppanggolin.geneFamily.GeneFamily`
        """
        return self._famGetter[name]

    def addGeneFamily(self, name):
        """
            Get the :class:`ppanggolin.geneFamily.GeneFamily` object that has the given `name`. If it does not exist,
            creates it.
            returns the geneFamily object.

            :param name: The gene family name to get if it exists, and create otherwise.
            :type name: str
        """
        fam = self._famGetter.get(name)
        if fam is None:
            fam = self._createGeneFamily(name)
        return fam

    """Graph methods"""
    @property
    def edges(self):
        """returns all the edges in the pangenome graph
        
        :return: list of :class:`ppanggolin.pangenome.Edge`
        :rtype: list
        """
        return list(self._edgeGetter.values())

    def addEdge(self, gene1, gene2):
        """
        Adds an edge between the two gene families that the two given genes belong to. Genes object are expected,
        and they are also expected to have a family assigned

        :param gene1: The first gene
        :type gene1: :class:`ppanggolin.genome.Gene`
        :param gene2: The second gene
        :type gene2: :class:`ppanggolin.genome.Gene`
        :return: the created Edge
        :rtype: :class:`ppanggolin.pangenome.Edge`
        """
        key = frozenset([gene1.family, gene2.family])
        edge = self._edgeGetter.get(key)
        if edge is None:
            edge = Edge(gene1, gene2)
            self._edgeGetter[key] = edge
        else:
            edge.addGenes(gene1, gene2)
        return edge

    """Organism methods"""
    @property
    def organisms(self):
        """returns all the organisms in the pangenome
        
        :return: list of :class:`ppanggolin.genome.Organism`
        :rtype: list
        """
        return list(self._orgGetter.values())

    def number_of_organisms(self):
        """Returns the number of organisms present in the pangenome
        
        :return: the number of organism
        :rtype: int
        """
        return len(self._orgGetter)

    def getOrganism(self, orgName):
        """
        Get an organism that is expected to be in the pangenome using its name, which is supposedly unique.
        Raises an error if the organism does not exist.

        :param orgName: Name of the :class:`ppanggolin.genome.Organism` to get
        :type orgName: str
        :return: The related Organism object
        :rtype: :class:`ppanggolin.genome.Organism`
        :raises KeyError: If the provided name is not in the pangenome
        """
        try:
            return self._orgGetter[orgName]
        except KeyError:
            raise KeyError(f"{orgName} does not seem to be in your pangenome")

    def addOrganism(self, newOrg):
        """
        adds an organism that did not exist previously in the pangenome if an :class:`ppanggolin.genome.Organism`
        object is provided. If an organism with the same name exists it will raise an error.
        If a :class:`str` object is provided, will return the corresponding organism that has this name
        OR create a new one if it does not exist.

        :param newOrg: Organism to add to the pangenome
        :type newOrg: :class:`ppanggolin.genome.Organism` or str
        :return: The created organism
        :rtype: :class:`ppanggolin.genome.Organism`
        :raises TypeError: if the provided `newOrg` is neither a str nor a :class:`ppanggolin.genome.Organism`
        """
        if isinstance(newOrg, Organism):
            oldLen = len(self._orgGetter)
            self._orgGetter[newOrg.name] = newOrg
            if len(self._orgGetter) == oldLen:
                raise KeyError(
                    f"Redondant organism name was found ({newOrg.name}). All of your organisms must have unique names.")
        elif isinstance(newOrg, str):
            org = self._orgGetter.get(newOrg)
            if org is None:
                org = Organism(newOrg)
                self._orgGetter[org.name] = org
            newOrg = org
        else:
            raise TypeError("Provide an Organism object or a str that will serve as organism name")
        return newOrg

    def getIndex(self):  # will not make a new index if it exists already
        """Creates an index for Organisms (each organism is assigned an Integer).

        :return: A dictionary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if not hasattr(self, "_orgIndex"):  # then the index does not exist yet
            self._orgIndex = {}
            for index, org in enumerate(self.organisms):
                self._orgIndex[org] = index
        return self._orgIndex

    def computeFamilyBitarrays(self, part='all'):
        """Based on the index generated by :meth:`ppanggolin.pangenome.Pangenome.getIndex`, generated a bitarray
        for each gene family.
        If the family j is present in the organism with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :return: A dictionnary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if not hasattr(self, "_orgIndex"):
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.getIndex()
        for fam in self.geneFamilies:
            fam.mkBitarray(self._orgIndex, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._orgIndex

    def get_fam_index(self):  # will not make a new index if it exists already
        """Creates an index for gene families (each family is assigned an Integer).

        :return: A dictionary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if not hasattr(self, "_fam_index"):  # then the index does not exist yet
            self._fam_index = {}
            for index, fam in enumerate(self.geneFamilies):
                self._fam_index[fam] = index
        return self._fam_index

    def compute_org_bitarrays(self, part='all'):
        """Based on the index generated by :meth:`ppanggolin.pangenome.Pangenome.get_fam_index`, generated a bitarray
        for each gene family.
        If the family j is present in the organism with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :return: A dictionary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if not hasattr(self, "_fam_index"):
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.get_fam_index()
        for org in self.organisms:
            org.mk_bitarray(index=self._fam_index, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._fam_index

    """RGP methods"""
    @property
    def regions(self):
        """returns all the regions (RGP) in the pangenome

        :return: list of :class:`ppanggolin.region.Region`
        :rtype: list
        """
        return list(self._regionGetter.values())

    def getOrAddRegion(self, regionName):
        """Returns a region with the given `regionName`. Creates it if it does not exist.

        :param regionName: The name of the region to return
        :type regionName: str
        :return: The region
        :rtype: :class:`ppanggolin.region.Region`
        """
        try:
            return self._regionGetter[regionName]
        except KeyError:  # then the region is not stored in this pangenome.
            newRegion = Region(regionName)
            self._regionGetter[regionName] = newRegion
            return newRegion

    def get_multigenics(self, dup_margin, persistent=True):
        """
        Returns the multigenic persistent families of the pangenome graph. A family will be considered multigenic
        if it is duplicated in more than `dup_margin` of the genomes where it is present.

        :param dup_margin: the ratio of presence in multicopy above which a gene family is considered multigenic
        :type dup_margin: float
        :param persistent: if we consider only the persistent genes
        :type persistent: bool
        :return: a `set` of gene families considered multigenic
        :rtype: set[:class:`ppanggolin.geneFamily.GeneFamily`]
        """
        multigenics = set()
        for fam in self.geneFamilies:
            if fam.namedPartition == "persistent" or not persistent:
                dup = len([genes for org, genes in fam.getOrgDict().items() if
                           len([gene for gene in genes if not gene.is_fragment]) > 1])
                if (dup / len(fam.organisms)) >= dup_margin:  # tot / nborgs >= 1.05
                    multigenics.add(fam)
        # logging.getLogger().info(f"{len(multigenics)} gene families are defined as being multigenic.
        # (duplicated in more than {dup_margin} of the genomes)")
        return multigenics

    def addRegions(self, regionGroup):
        """Takes an Iterable or a Region object and adds it to the pangenome

        :param regionGroup: a region or an Iterable of regions to add to the pangenome
        :type regionGroup: :class:`ppanggolin.region.Region` or Iterable[:class:`ppanggolin.region.Region`]
        :raises TypeError: if regionGroup is neither a Region nor a Iterable[:class:`ppanggolin.region.Region`]
        """
        oldLen = len(self._regionGetter)
        if isinstance(regionGroup, Iterable):
            for region in regionGroup:
                self._regionGetter[region.name] = region
            if len(self._regionGetter) != len(regionGroup) + oldLen:
                raise Exception("Two regions had an identical name, which was unexpected.")
        elif isinstance(regionGroup, Region):
            self._regionGetter[regionGroup.name] = regionGroup
        else:
            raise TypeError(f"An iterable or a 'Region' type object were expected, "
                            f"but you provided a {type(regionGroup)} type object")

    """Spot methods"""
    def addSpots(self, spots):
        """Adds the given iterable of spots to the pangenome.

        :param spots: An iterable of :class:`ppanggolin.region.Spot`.
        :type spots: Iterable[:class:`ppanggolin.region.Spot`]
        """
        self.spots |= set(spots)

    """Modules methods"""
    def addModules(self, modules):
        """Adds the given iterable of modules to the pangenome

        :param modules: an iterable of :class:`ppanggolin.module.Module`
        :type modules: Iterable[:class:`ppanggolin.module.Module`]
        """
        self.modules |= set(modules)

    def compute_mod_bitarrays(self, part='all'):
        """Based on the index generated by :meth:`ppanggolin.pangenome.Pangenome.get_fam_index`, generated a bitarray
        for each gene family.
        If the family j is present in the module with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :return: A dictionary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if not hasattr(self, "_fam_index"):
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.get_fam_index()
        for module in self.modules:
            module.mk_bitarray(index=self._fam_index, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._fam_index
