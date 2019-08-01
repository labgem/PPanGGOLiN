#!/usr/bin/env python3
#coding: utf8

from collections import defaultdict

class RNA:
    def __init__(self, ID):
        self.ID = ID
        self.is_fragment = False
        self.type = ""

    def fill_annotations(self, start, stop, strand, geneType = "", name = "", product=""):
        self.start = int(start)
        self.stop = int(stop)
        self.type = geneType
        self.strand = strand
        self.product = product
        self.name = name

    def fill_parents(self, organism, contig):
        self.organism = organism
        self.contig = contig

    def add_dna(self, dna):
        if not isinstance(dna, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(dna)}' type object")
        self.dna = dna


class Gene:
    def __init__(self, ID):
        self.ID = ID
        self.is_fragment = False
        self.type = ""
        self.position = None
        self.family = None

    def fill_annotations(self, start, stop, strand, geneType = "", position = None, name = "", product="", genetic_code = 11):
        self.start = int(start)
        self.stop = int(stop)
        self.type = geneType
        self.strand = strand
        self.position = position
        self.product = product
        self.name = name
        self.genetic_code = genetic_code

    def fill_parents(self, organism, contig):
        #check the instances ???
        self.organism = organism
        self.contig = contig

    def add_protein(self, protein):
        if not isinstance(protein, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(protein)}' type object")
        self.protein = protein

    def add_dna(self, dna):
        if not isinstance(dna, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(dna)}' type object")
        self.dna = dna

class Contig:
    def __init__(self, name, is_circular = False):
        self.name = name
        self.is_circular = is_circular
        self.RNAs = set()#saving the rna annotations. We're not using them in the vast majority of cases.
        self._genes_start = {}
        self._genes_position = []
    
    @property
    def genes(self):
        return self._genes_start.values()

    def __str__(self):
        return self.name

    def addSequence(self, seq):
        if not isinstance(seq, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(seq)}' type object")
        self.sequence = seq

    def __iter__(self):
        return iter(self.genes)

    # retrieve gene by start position
    def __getitem__(self, index):
        gene = self._genes_start.get(index)
        if not gene:
            if type(index) != int:
                raise TypeError(f"Expected type is int, given type was '{type(index)}'")
            raise IndexError(f"No gene start at the given position {index}")
        return gene
    
    def addRNA(self, gene):
        if not isinstance(gene, RNA):
            raise TypeError(f"'Gene' type was expected but you provided a '{type(gene)}' type object")
        self.RNAs.add(gene)

    def addGene(self, gene):
        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type was expected but you provided a '{type(gene)}' type object")
        if gene.position is None:
            raise TypeError(f"The gene object needs to have its position in the contig filled before adding it")
        while len(self._genes_position) <= gene.position:
            # adding empty values. They should be filled by the end of the parsing. Doing this because genes are not always met in order.
            self._genes_position.append(None)
        self._genes_position[gene.position] = gene
        self._genes_start[gene.start] = gene

class Organism:
    def __init__(self, name):
        self.name = name
        self._contigs_getter = {}

    @property
    def contigs(self):
        return set(self._contigs_getter.values())

    def __iter__(self):
        return iter([ gene for gene in contig for contig in self._contigs_getter.values()])

    def __str__(self):
        return self.name

    def copy(self):
        """ returns a new instance with the same attributes and values than self, identical contigs, but some changes in the genes"""
        new_org = Organism(self.name)
        for contigObj in self.contigs:
            new_contig = Contig(contigObj.name, is_circular=contigObj.is_circular)
            new_org._contigs_getter[new_contig.name] = new_contig
            for gene in contigObj.genes:
                new_gene= Gene(gene.ID)
                new_gene.fill_annotations(gene.start, gene.stop, gene.strand, gene.type, gene.position, gene.name, gene.product, gene.genetic_code)
                new_gene.fill_parents(new_org, new_contig)
                new_gene.family = gene.family
                new_gene.is_fragment = gene.is_fragment
                new_contig.addGene(new_gene)
        return new_org

    def addContig(self, key, is_circular = False):
        contig = self._contigs_getter.get(key)
        if contig is None:
            contig = self._createContig(key, is_circular)
        return contig

    def _createContig(self, key, is_circular = False):
        new_contig = Contig(key, is_circular)
        self._contigs_getter[key] = new_contig
        return new_contig