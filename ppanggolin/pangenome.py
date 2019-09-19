#!/usr/bin/env python3
#coding: utf8

#default libraries
from collections import defaultdict
#installed libraries
import gmpy2

#local libraries
from ppanggolin.genome import Organism, Gene

class Region:
    def __init__(self, ID):
        self.genes = []
        self.name = ID
        self.score = 0

    def append(self, value):
        # allowing only gene-class objects in a region.
        if isinstance(value, Gene):
            self.genes.append(value)
        else:
            raise TypeError("Unexpected class / type for " + type(value) +
                            " when adding it to a region of genomic plasticity")
    @property
    def organism(self):
        return self.genes[0].organism

    @property
    def contig(self):
        return self.genes[0].contig

    def __len__(self):
        return len(self.genes)

    def __getitem__(self, index):
        return self.genes[index]

    def __iter__(self):
        return iter(self.genes)

    def __str__(self):
        return self.name

class Edge:
    def __init__(self, sourceGene, targetGene):
        if sourceGene.family is None:
            raise Exception(f"You cannot create a graph without gene families. gene {sourceGene.ID} did not have a gene family.")
        if targetGene.family is None:
            raise Exception(f"You cannot create a graph without gene families. gene {targetGene.ID} did not have a gene family.")
        self.source = sourceGene.family
        self.target = targetGene.family
        self.source._edges[self.target] = self
        self.target._edges[self.source] = self
        self.organisms = defaultdict(list)
        self.addGenes(sourceGene, targetGene)

    def getOrgDict(self):
        return self.organisms

    @property
    def genePairs(self):
        return [ gene_pair for gene_list in self.organisms.values() for gene_pair in gene_list ]

    def addGenes(self, sourceGene, targetGene):
        org = sourceGene.organism
        if org != targetGene.organism:
            raise Exception(f"You tried to create an edge between two genes that are not even in the same organism ! (genes are '{sourceGene.ID}' and '{targetGene.ID}')")
        self.organisms[org].append((sourceGene, targetGene))

    def remove(self):
        self.source._edges[self.target].discard(self)
        self.target._edges[self.source].discard(self)

class GeneFamily:
    def __init__(self, ID, name):
        self.name = name
        self.ID = ID
        self._edges = {}
        self._genePerOrg = defaultdict(set)
        self.genes = set()
        self.removed = False#for the repeated family not added in the main graph
        self.sequence = ""
        self.partition = ""

    def addSequence(self, seq):
        self.sequence = seq

    def addPartition(self, partition):
        self.partition = partition

    @property
    def namedPartition(self):
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

    def addGene(self, gene):
        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type object was expected, but '{type(gene)}' type object was provided.")
        self.genes.add(gene)
        gene.family = self
        if hasattr(gene, "organism"):
            self._genePerOrg[gene.organism].add(gene)

    def mkBitarray(self, index):
        """ produces a bitarray representing the presence / absence of the family in the pangenome"""
        self.bitarray = gmpy2.xmpz(0)
        for org in self.organisms:
            self.bitarray[index[org]] = 1

    def getOrgDict(self):
        try:
            return self._genePerOrg
        except AttributeError:
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg

    def getGenesPerOrg(self, org):
        try:
            return self._genePerOrg[org]
        except AttributeError:
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg[org]

    @property
    def neighbors(self):
        return set(self._edges.keys())

    @property
    def edges(self):
        return self._edges.values()

    @property
    def organisms(self):
        try:
            return self._genePerOrg.keys()
        except AttributeError:#then the genes have been added before they had organisms
            for gene in self.genes:
                self._genePerOrg[gene.organism].add(gene)
            return self._genePerOrg.keys()

class Pangenome:
    def __init__(self):
        #basic parameters
        self._famGetter = {}
        self.max_fam_id = 0
        self._orgGetter = {}
        self._edgeGetter = {}
        self.status = {
                    'genomesAnnotated': "No",
                    'geneSequences' : "No",
                    'genesClustered':  "No",
                    'defragmented':"No",
                    'geneFamilySequences':"No",
                    'neighborsGraph':  "No",
                    'partitionned':  "No",
                    'predictedRGP' : "No"
                }
        self.parameters = {}

    def addFile(self, pangenomeFile):
        from ppanggolin.formats import getStatus#importing on call instead of importing on top to avoid cross-reference problems.
        getStatus(self, pangenomeFile)
        self.file = pangenomeFile

    @property
    def genes(self):
        if len(self.organisms) > 0:#if we have organisms, they're supposed to have genes
            return [ gene for org in self._orgGetter.values() for contig in org.contigs for gene in contig.genes ]
        elif len(self.geneFamilies) > 0:#else, the genes will be stored in the gene families (maybe)
            return [ gene for geneFam in self.geneFamilies for gene in geneFam.genes ]

    @property
    def geneFamilies(self):
        return self._famGetter.values()

    @property
    def edges(self):
        return self._edgeGetter.values()

    @property
    def organisms(self):
        return self._orgGetter.values()

    def number_of_organisms(self):
        return len(self._orgGetter)

    def number_of_geneFamilies(self):
        return len(self._famGetter)

    def _yield_genes(self):
        """
            Use a generator to get all the genes of a pangenome
        """
        if self.number_of_organisms() > 0:#if we have organisms, they're supposed to have genes
            for org in self.organisms:
                for contig in org.contigs:
                     for gene in contig.genes:
                         yield gene
        elif self.number_of_geneFamilies() > 0:
            for geneFam in self.geneFamilies:
                for gene in geneFam.genes:
                    yield gene

    def _mkgeneGetter(self):
        """
            Since the genes are never explicitely 'added' to a pangenome (but rather to a gene family, or a contig), the pangenome cannot directly extract a gene from a geneID since it does not 'know' them.
            if at some point we want to extract genes from a pangenome we'll create a geneGetter.
            The assumption behind this is that the pangenome has been filled and no more gene will be added.
        """
        self._geneGetter = {}
        for gene in self._yield_genes():
            self._geneGetter[gene.ID] = gene

    def getGene(self, geneID):
        try:
            return self._geneGetter[geneID]
        except AttributeError:#in that case, either the gene getter has not been computed, or the geneID is not in the pangenome.
            self._mkgeneGetter()#make it
            return self.getGene(geneID)#return what was expected. If the geneID does not exist it will raise an error.
        except KeyError:
            return None

    def info(self):
        infostr = ""
        infostr += f"Gene families : {len(self.geneFamilies)}\n"
        infostr += f"Organisms : {len(self.organisms)}\n"
        nbContig = 0
        for org in self.organisms:
            for _ in org.contigs:
                nbContig+=1
        infostr += f"Contigs : {nbContig}\n"
        infostr += f"Genes : {len(self.genes)}\n"
        infostr += f"Edges : {len(self.edges)}\n"
        nbP=0
        nbC=0
        nbS=0
        for fam in self.geneFamilies:
            if fam.partition == "C":
                nbC+=1
            elif fam.partition == "P":
                nbP+=1
            elif fam.partition.startswith("S"):
                nbS+=1
        infostr += f"Persistent : {nbP}\n"
        infostr += f"Shell : {nbS}\n"
        infostr += f"Cloud : {nbC}\n"

        return infostr

    def addOrganism(self, newOrg):
        """
            adds an organism that did not exist previously in the pangenome if an Organism object is provided.
            If a str object is provided, will return the corresponding organism OR create a new one.
        """
        if isinstance(newOrg, Organism):
            oldLen = len(self._orgGetter)
            self._orgGetter[newOrg.name] = newOrg
            if len(self._orgGetter) == oldLen:
                raise KeyError(f"Redondant organism name was found ({newOrg.name}). All of your organisms must have unique names.")
        elif isinstance(newOrg, str):
            org = self._orgGetter.get(newOrg)
            if org is None:
                org = Organism(newOrg)
                self._orgGetter[newOrg.name] = newOrg
            newOrg = org
        return newOrg

    def addGeneFamily(self, name):
        """
            Creates a geneFamily object with the provided name and adds it to the pangenome if it does not exist.
            Otherwise, does not create anything.
            returns the geneFamily object.
        """
        fam = self._famGetter.get(name)
        if fam is None:
            fam = self._createGeneFamily(name)
        return fam

    def getGeneFamily(self, name):
        return self._famGetter[name]

    def addEdge(self, gene1, gene2):
        key = frozenset([gene1.family,gene2.family])
        edge = self._edgeGetter.get(key)
        if edge is None:
            edge = Edge(gene1, gene2)
            self._edgeGetter[key] = edge
        else:
            edge.addGenes(gene1,gene2)
        return edge

    def _createGeneFamily(self, name):
        newFam = GeneFamily(ID = self.max_fam_id, name = name)
        self.max_fam_id+=1
        self._famGetter[newFam.name] = newFam
        return newFam

    def getIndex(self):#will not make a new index if it exists already
        if not hasattr(self, "_orgIndex"):#then the index does not exist yet
            self._orgIndex = {}
            for index, org in enumerate(self.organisms):
                self._orgIndex[org] = index
        return self._orgIndex

    def computeFamilyBitarrays(self):
        if not hasattr(self, "_orgIndex"):#then the bitarrays don't exist yet, since the org index does not exist either.
            self.getIndex()
            for fam in self.geneFamilies:
                fam.mkBitarray(self._orgIndex)
        #case where there is an index but the bitarrays have not been computed???
        return self._orgIndex