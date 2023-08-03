#!/usr/bin/env python3
# coding: utf8

# default libraries
from collections import defaultdict
from typing import Dict, Generator, List, Tuple

from ppanggolin.genome import Gene, Organism


class Edge:
    """The Edge class represents an edge between two gene families in the pangenome graph. It is associated with all the
     organisms in which the neighborship is found, and all the involved genes as well.
    Methods:
    - __init__(self, source_gene: Gene, target_gene: Gene): Constructor method that initializes an Edge object with a source gene and a target gene.
    - get_org_dict(self) -> Dict[Organism, List[Tuple[Gene, Gene]]]: Returns a dictionary with organisms as keys and an iterable of the pairs of genes as values.
    - gene_pairs(self) -> List[Tuple[Gene, Gene]]: Returns a list of all the gene pairs of the Edge.
    - add_genes(self, source_gene: Gene, target_gene: Gene): Adds genes to the edge. They are supposed to be on the same organism.

    Fields:
    - source: A GeneFamily object representing the source gene family of the edge.
    - target: A GeneFamily object representing the target gene family of the edge.
    - organisms: A defaultdict object representing the organisms in which the edge is found and the pairs of genes involved.
    """

    def __init__(self, source_gene: Gene, target_gene: Gene):
        """Constructor method

        :param source_gene: a first gene to initialize the edge
        :param target_gene: a second gene to initialize the edge
        """
        # TODO try to change for gene family ?
        if source_gene.family is None:
            raise AttributeError(f"You cannot create a graph without gene families. "
                                 f"gene {source_gene.ID} did not have a gene family.")
        if target_gene.family is None:
            raise AttributeError(f"You cannot create a graph without gene families. "
                                 f"gene {target_gene.ID} did not have a gene family.")
        self.source = source_gene.family
        self.target = target_gene.family
        self.source.set_edge(self.target, self)
        self.target.set_edge(self.source, self)
        self._organisms = defaultdict(list)
        self.add_genes(source_gene, target_gene)

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """Get all the organisms belonging to the edge

        :return: Generator with organisms as key and an iterable of the pairs of genes as value
        """
        for organism in self._organisms.keys():
            yield organism

    @property
    def number_of_organisms(self):
        return len(self._organisms)

    def get_organism_genes_pairs(self, organism: Organism):
        return self._organisms[organism]

    def get_organisms_dict(self):
        return self._organisms

    @property
    def gene_pairs(self) -> List[Tuple[Gene, Gene]]:
        """ Get list of all the gene pairs of the Edge

        :return: A list of all the gene pairs of the Edge
        """
        return [gene_pair for gene_list in self.get_organisms_dict().values() for gene_pair in gene_list]

    def add_genes(self, source_gene: Gene, target_gene: Gene):
        """Adds genes to the edge. They are supposed to be on the same organism.

        :param source_gene: a source gene to add to the edge
        :param target_gene: a target gene to add to the edge

        :raises TypeError: If the genes are not with Gene type
        :raises ValueError: If genes are not associated to an organism
        :raises Exception: If the genes are not on the same organism.
        """
        if not isinstance(source_gene, Gene) or not isinstance(target_gene, Gene):
            raise TypeError(f"Genes are expected to be added to edge. "
                            f"Given type for source: {type(source_gene)} and target: {type(target_gene)}")
        if source_gene.organism is None or target_gene.organism is None:
            raise ValueError("Genes are not associated to organism. It's needed to create add genes to edge")
        if source_gene.organism != target_gene.organism:
            raise Exception(f"You tried to create an edge between two genes that are not even in the same organism ! "
                            f"(genes are '{source_gene.ID}' and '{target_gene.ID}')")
        self._organisms[source_gene.organism].append((source_gene, target_gene))
