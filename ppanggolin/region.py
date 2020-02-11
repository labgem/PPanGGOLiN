#!/usr/bin/env python3
#coding: utf8

#default libraries
import logging
from collections.abc import Iterable

#local libraries
from ppanggolin.genome import Organism, Gene

class Region:
    def __init__(self, ID):
        self.genes = []
        self.name = ID
        self.score = 0

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        """ expects another Region type object. Will test whether two Region objects have the same gene families"""
        if not isinstance(other, Region):
            raise TypeError(f"'Region' type object was expected, but '{type(other)}' type object was provided.")
        if [ gene.family for gene in self.genes ] == [ gene.family for gene in other.genes ]:
            return True
        if [ gene.family for gene in self.genes ] == [ gene.family for gene in other.genes[::-1]]:
            return True
        return False

    def append(self, value):
        # allowing only gene-class objects in a region.
        if isinstance(value, Gene):
            self.genes.append(value)
        else:
            raise TypeError("Unexpected class / type for " + type(value) +" when adding it to a region of genomic plasticity")

    @property
    def families(self):
        return { gene.family for gene in self.genes }

    @property
    def start(self):
        return min(self.genes, key = lambda x : x.start).start

    @property
    def startGene(self):
        return min(self.genes, key = lambda x : x.position)

    @property
    def stopGene(self):
        return max(self.genes, key = lambda x : x.position)

    @property
    def stop(self):
        return max(self.genes, key = lambda x : x.stop).stop

    @property
    def organism(self):
        return self.genes[0].organism

    @property
    def contig(self):
        return self.genes[0].contig

    @property
    def isWholeContig(self):
        """ Indicates if the region is an entire contig """
        if self.startGene.position == 0 and self.stopGene.position == len(self.contig.genes)-1:
            return True
        return False

    @property
    def isContigBorder(self):
        if len(self.genes) == 0:
            raise Exception("Your region has no genes. Something wrong happenned.")
        if self.startGene.position == 0 and not self.contig.is_circular:
            return True
        elif self.stopGene.position == len(self.contig.genes)-1 and not self.contig.is_circular:
            return True
        return False

    def getRNAs(self):
        RNAs = set()
        for rna in self.contig.RNAs:
            if rna.start > self.start and rna.start < self.stop:
                RNAs.add(rna)
        return RNAs

    def __len__(self):
        return len(self.genes)

    def __getitem__(self, index):
        return self.genes[index]

    def getBorderingGenes(self, n, multigenics):
        border = [[], []]
        pos = self.startGene.position
        init = pos
        while len(border[0]) < n and (pos != 0 and not self.contig.is_circular):
            curr_gene = None
            if pos == 0:
                if self.contig.is_circular:
                    curr_gene = self.contig.genes[-1]
            else:
                curr_gene = self.contig.genes[pos -1]
            if curr_gene is not None and curr_gene.family not in multigenics and curr_gene.family.namedPartition == "persistent":
                border[0].append(curr_gene)
            pos -= 1
            if pos == -1 and self.contig.is_circular:
                pos = len(self.contig.genes)
            if pos == init:
                logging.getLogger().warning("looped around the contig")
                break#looped around the contig
        pos = self.stopGene.position
        init = pos
        while len(border[1]) < n and (pos != len(self.contig.genes)-1 and not self.contig.is_circular):
            curr_gene = None
            if pos == len(self.contig.genes)-1:
                if self.contig.is_circular:
                    curr_gene = self.contig.genes[0]
            else:
                curr_gene = self.contig.genes[pos+1]
            if curr_gene is not None and curr_gene.family not in multigenics:
                border[1].append(curr_gene)
            pos+=1
            if pos == len(self.contig.genes) and self.contig.is_circular:
                pos = -1
            if pos == init:
                logging.getLogger().warning("looped around the contig")
                break#looped around the contig
        return border


class Spot:
    def __init__(self, ID):
        self.ID = ID
        self.regions = set()
        self._uniqOrderedSet = {}
        self._compOrderedSet = False
        self._uniqContent = {}
        self._compContent = False

    def addRegions(self, regions):
        """ Adds region(s) contained in an Iterable to the spot which all have the same bordering persistent genes provided with 'borders'"""
        if isinstance(regions, Iterable):
            for region in regions:
                self.addRegion(region)
        else:
            raise Exception("The provided 'regions' variable was not an Iterable")

    def addRegion(self, region):
        if isinstance(region, Region):
            self.regions.add(region)

    def borders(self, set_size, multigenics):
        """ extracts all the borders of all RGPs belonging to the spot"""
        raise NotImplementedError()

    def _mkUniqOrderedSetObj(self):
        """cluster RGP into groups that have an identical synteny""" 
        for rgp in self.regions:
            z = True
            for seenRgp in self._uniqOrderedSet:
                if rgp == seenRgp:
                    z = False
                    self._uniqOrderedSet[seenRgp].add(rgp)
            if z:
                self._uniqOrderedSet[rgp] = set([rgp])

    def _mkUniqContent(self):
        """cluster RGP into groups that have identical gene content"""
        for rgp in self.regions:
            z = True
            for seenRgp in self._uniqContent:
                if rgp.families == seenRgp.families:
                    z = False
                    self._uniqContent[seenRgp].add(rgp)
            if z:
                self._uniqContent[rgp] = set([rgp])

    def _getContent(self):
        """Creates the _uniqContent object if it was never computed. Return it in any case"""
        if not self._compContent:
            self._mkUniqContent()
            self._compContent = True
        return self._uniqContent

    def _getOrderedSet(self):
        """Creates the _uniqSyn object if it was never computed. Return it in any case"""
        if not self._compOrderedSet:
            self._mkUniqOrderedSetObj()
            self._compOrderedSet = True
        return self._uniqOrderedSet

    def getUniqOrderedSet(self):
        """ returns an Iterable of all the unique syntenies in the spot"""
        return set(self._getOrderedSet().keys())

    def getUniqContent(self):
        """ returns an Iterable of all the unique rgp (in terms of gene family content) in the spot"""
        return set(self._getContent().keys())

    def countUniqContent(self):
        """Returns a counter with a representative rgp as key and the number of identical rgp in terms of gene family content as value"""
        return dict([ (key, len(val)) for key, val in self._getContent().items()])

    def countUniqOrderedSet(self):
        """ Returns a counter with a representative rgp as key and the number of identical rgp in terms of synteny as value"""
        return dict([ (key, len(val)) for key, val in self._getOrderedSet().items()])