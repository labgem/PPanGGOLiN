#!/usr/bin/env python3

# default libraries
from collections import defaultdict
from typing import Dict, Generator, List, Tuple

from ppanggolin.genome import Gene, Organism


class Edge:
    """
    The Edge class represents an edge between two gene families in the pangenome graph. It is associated with all the
    organisms in which the neighborship is found, and all the involved genes as well.

    Methods:
        - get_org_dict: Returns a dictionary with organisms as keys and an iterable of the pairs in genes as values.
        - gene_pairs: Returns a list of all the gene pairs in the Edge.
        - add_genes: Adds genes to the edge. They are supposed to be in the same organism.

    Fields:
        - source: A GeneFamily object representing the source gene family of the edge.
        - target: A GeneFamily object representing the target gene family of the edge.
        - organisms: A defaultdict object representing the organisms in which the edge is found and the pairs of genes involved.
    """

    def __init__(self, source_gene: Gene, target_gene: Gene):
        """Constructor method

        :param source_gene: First gene to initialize the edge
        :param target_gene: Second gene to initialize the edge
        """
        # TODO try to change for gene family ?
        if source_gene.family is None:
            raise AttributeError(
                f"You cannot create a graph without gene families. "
                f"gene {source_gene.ID} did not have a gene family."
            )
        if target_gene.family is None:
            raise AttributeError(
                f"You cannot create a graph without gene families. "
                f"gene {target_gene.ID} did not have a gene family."
            )
        self.source = source_gene.family
        self.target = target_gene.family
        self.source.set_edge(self.target, self)
        self.target.set_edge(self.source, self)
        self._organisms = defaultdict(list)
        self.add_genes(source_gene, target_gene)

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """Get all the organisms belonging to the edge

        :return: Generator with organisms as the key and an iterable of the gene pairs as value
        """
        yield from self._organisms.keys()

    @property
    def number_of_organisms(self) -> int:
        """Get the number of organisms in the edge

        :return: Number of organisms
        """
        return len(self._organisms)

    def get_organism_genes_pairs(self, organism: Organism) -> List[Tuple[Gene, Gene]]:
        """Get the gene pair corresponding to the given organism

        :param organism: Wanted organism

        :return: Pair of genes in the edge corresponding to the given organism
        """
        return self._organisms[organism]

    def get_organisms_dict(self) -> Dict[Organism, List[Tuple[Gene, Gene]]]:
        """Get all the organisms with their corresponding pair of genes in the edge

        :return: Dictionary with the organism as the key and list of gene pairs as value
        """
        return self._organisms

    @property
    def gene_pairs(self) -> List[Tuple[Gene, Gene]]:
        """Get the list of all the gene pairs in the Edge

        :return: A list of all the gene pairs in the Edge
        """
        return [
            gene_pair
            for gene_list in self.get_organisms_dict().values()
            for gene_pair in gene_list
        ]

    def add_genes(self, source_gene: Gene, target_gene: Gene):
        """
        Adds genes to the edge.
        They are supposed to be in the same organism.

        :param source_gene: Gene corresponding to the source of the edge
        :param target_gene: Gene corresponding to the target of the edge

        :raises TypeError: If the genes are not with Gene type
        :raises ValueError: If genes are not associated with an organism
        :raises Exception: If the genes are not in the same organism.
        """
        if not isinstance(source_gene, Gene) or not isinstance(target_gene, Gene):
            raise TypeError(
                f"Genes are expected to be added to edge. "
                f"Given type for source: {type(source_gene)} and target: {type(target_gene)}"
            )
        if source_gene.organism is None or target_gene.organism is None:
            raise ValueError(
                "Genes are not associated to genome. It's needed to create add genes to edge"
            )
        if source_gene.organism != target_gene.organism:
            raise Exception(
                f"You tried to create an edge between two genes that are not even in the same genome ! "
                f"(genes are '{source_gene.ID}' and '{target_gene.ID}')"
            )
        self._organisms[source_gene.organism].append((source_gene, target_gene))
