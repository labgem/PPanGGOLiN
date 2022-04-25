#!/usr/bin/env python3
# coding: utf8

# installed libraries
import logging

import gmpy2


class Feature:
    def __init__(self, identifier):
        self.ID = identifier
        self.is_fragment = False
        self.type = ""
        self.start = None
        self.stop = None
        self.type = None
        self.strand = None
        self.product = None
        self.name = None
        self.local_identifier = None
        self.organism = None
        self.contig = None
        self.dna = None

    def fill_annotations(self, start, stop, strand, gene_type="", name="", product="", local_identifier="",
                         position=None, genetic_code=11):
        # genetic code, and position are not used in the default function.
        self.start = int(start)
        self.stop = int(stop)
        self.type = gene_type
        self.strand = strand
        self.product = product
        self.name = name
        self.local_identifier = local_identifier

    def fill_parents(self, organism, contig):
        self.organism = organism
        self.contig = contig

    def add_dna(self, dna):
        if not isinstance(dna, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(dna)}' type object")
        self.dna = dna


class RNA(Feature):
    pass


class Gene(Feature):
    def __init__(self, gene_id):
        super().__init__(gene_id)
        self.position = None
        self.family = None
        self.RGP = set()
        self.genetic_code = None
        self.protein = None

    def __str__(self):
        return str(self.ID)

    def fill_annotations(self, start, stop, strand, gene_type="", name="", product="", local_identifier="",
                         position=None, genetic_code=11):
        super().fill_annotations(start, stop, strand, gene_type, name, product, local_identifier)
        self.position = position
        self.genetic_code = genetic_code

    def add_protein(self, protein):
        if not isinstance(protein, str):
            raise TypeError(f"'str' type was expected but you provided a '{type(protein)}' type object")
        self.protein = protein


class Contig:
    def __init__(self, name, is_circular=False):
        self.name = name
        self.is_circular = is_circular
        self.RNAs = set()  # saving the rna annotations. We're not using them in the vast majority of cases.
        self._genes_start = {}
        self._genes_position = []

    @property
    def genes(self):
        return self._genes_position

    def __str__(self):
        return self.name

    def __iter__(self):
        return iter(self.genes)

    # retrieve gene by start position
    def __getitem__(self, index):
        gene = self._genes_start.get(index)
        if not gene:
            if not isinstance(index, int):
                raise TypeError(f"Expected type is int, given type was '{type(index)}'")
            raise IndexError(f"No gene start at the given position {index}")
        return gene

    def add_rna(self, gene):
        if not isinstance(gene, RNA):
            raise TypeError(f"'RNA' type was expected but you provided a '{type(gene)}' type object")
        self.RNAs.add(gene)

    def add_gene(self, gene):
        if not isinstance(gene, Gene):
            raise TypeError(f"'Gene' type was expected but you provided a '{type(gene)}' type object")
        if gene.position is None:
            raise TypeError(f"The gene object needs to have its position in the contig filled before adding it")
        while len(self._genes_position) <= gene.position:
            # adding empty values. They should be filled by the end of the parsing.
            # Doing this because genes are not always met in order.
            self._genes_position.append(None)
        self._genes_position[gene.position] = gene
        self._genes_start[gene.start] = gene


class Organism:
    def __init__(self, name):
        self.name = name
        self._contigs_getter = {}
        self.bitarray = None

    @property
    def families(self):
        """returns the gene families present in the organism"""
        return {gene.family for contig in self.contigs for gene in contig.genes}

    @property
    def genes(self):
        for contig in self.contigs:
            for gene in contig.genes:
                yield gene

    def number_of_genes(self):
        return sum([len(list(contig.genes)) for contig in self.contigs])

    @property
    def contigs(self):
        return self._contigs_getter.values()

    def __str__(self):
        return self.name

    def get_or_add_contig(self, key, is_circular=False):
        contig = self._contigs_getter.get(key)
        if contig is None:
            contig = self._create_contig(key, is_circular)
        return contig

    def _create_contig(self, key, is_circular=False):
        new_contig = Contig(key, is_circular)
        self._contigs_getter[key] = new_contig
        return new_contig

    def mk_bitarray(self, index, partition='all'):
        """Produces a bitarray representing the presence / absence of families in the organism using the provided index
        The bitarray is stored in the :attr:`bitarray` attribute and is a :class:`gmpy2.xmpz` type.
        :param partition: Filter partition
        type partition: str
        :param index: The index computed by :func:`ppanggolin.pangenome.Pangenome.getIndex`
        :type index: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        self.bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        if partition == 'all':
            logging.getLogger().debug(f"all")
            for fam in self.families:
                self.bitarray[index[fam]] = 1
        elif partition in ['shell', 'cloud']:
            logging.getLogger().debug(f"shell, cloud")
            for fam in self.families:
                if fam.named_partition == partition:
                    self.bitarray[index[fam]] = 1
        elif partition == 'accessory':
            logging.getLogger().debug(f"accessory")
            for fam in self.families:
                if fam.named_partition in ['shell', 'cloud']:
                    self.bitarray[index[fam]] = 1
        else:
            raise Exception("There is not any partition corresponding please report a github issue")
