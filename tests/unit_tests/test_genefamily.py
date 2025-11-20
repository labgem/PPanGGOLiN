#! /usr/bin/env python3

import pytest
from random import randint
from typing import Generator, Set
from itertools import combinations_with_replacement

from ppanggolin.pangenome import Edge
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.genome import Gene, Organism, Contig
from ppanggolin.region import Spot, Module


class TestGeneFamily:
    """Tests the gene family class"""

    @pytest.fixture
    def family(self) -> Generator[GeneFamily, None, None]:
        """Create a gene family for all tests"""
        yield GeneFamily(1, "test")

    def test_construct_gene_family(self, family):
        """Tests that a GeneFamily object can be created with valid family_id and name"""
        assert isinstance(family, GeneFamily)
        assert all(
            attr
            in [
                "ID",
                "name",
                "_edges_getter",
                "_genePerOrg",
                "_genes_getter",
                "_representative",
                "removed",
                "sequence",
                "_partition",
                "_spots",
                "_module",
                "bitarray",
                "_metadata_getter",
            ]
            for attr in family.__dict__
        )  # Check that no attribute was added else it should be tested
        assert all(
            hasattr(family, attr)
            for attr in [
                "ID",
                "name",
                "_edges_getter",
                "_genePerOrg",
                "_genes_getter",
                "removed",
                "sequence",
                "_partition",
                "_spots",
                "_module",
                "bitarray",
            ]
        )  # Check that no attribute was removed else it should be tested
        assert family.ID == 1
        assert family.name == "test"
        assert family._edges_getter == {}
        assert family._genePerOrg == {}
        assert family._genes_getter == dict()
        assert not family.removed  # for the repeated family not added in the main graph
        assert family.sequence == ""
        assert family.partition == ""
        assert family._spots == set()
        assert family._module is None
        assert family.bitarray is None

    @pytest.mark.parametrize(
        "partition, name",
        [
            ("P", "persistent"),
            ("Pp", "persistent"),
            ("P whatever, only first letter is important", "persistent"),
            ("C", "cloud"),
            ("C loud", "cloud"),
            ("C whatever, only first letter is important", "cloud"),
            ("S", "shell"),
            ("Shut", "shell"),
            ("S whatever, only first letter is important", "shell"),
            ("un de troa kvar", "undefined"),
            ("1", "undefined"),
            ("p", "undefined"),
            ("c", "undefined"),
            ("s", "undefined"),
        ],
    )
    def test_get_named_partition_of_gene_family_object(self, family, partition, name):
        """Tests that the named partition of a GeneFamily object can be retrieved"""
        family.partition = partition
        assert family.named_partition == name

    def test_get_named_partition_error_partition_empty(self, family):
        """Tests that if no partition given to gene family, raise a ValueError"""
        with pytest.raises(ValueError):
            _ = family.named_partition

    def test_add_sequence_to_gene_family(self, family):
        """Tests that a sequence can be added to a GeneFamily object"""
        family.add_sequence("ATCG")
        assert family.sequence == "ATCG"

    def test_add_gene_to_gene_family(self, family):
        """Tests that a Gene object can be added to a GeneFamily object"""
        gene = Gene("gene1")
        family.add(gene)
        assert gene in family.genes
        assert gene.family == family

    def test_add_gene_error(self, family):
        """Tests that a non-gene object can't be added to a GeneFamily as gene"""
        with pytest.raises(TypeError):
            family.add(33)

    def test_set_representative_gene(self, family):
        gene = Gene("representative_gene")
        family.representative = gene
        assert family._representative == gene

    def test_get_representative_gene(self, family):
        gene = Gene("representative_gene")
        family.representative = gene
        assert family.representative == gene

    def test_raise_typeerror_with_no_gene_type_as_representative(self, family):
        with pytest.raises(TypeError):
            family.representative = "test"

    def test_raise_exception_if_representative_not_set(self, family):
        with pytest.raises(Exception):
            _ = family.representative

    @pytest.fixture
    def genes(self) -> Generator[Set[Gene], None, None]:
        """Create a set of genes to fill gene families"""
        genes = set()
        for i in range(1, randint(11, 20)):
            gene = Gene(f"gene_{str(i)}")
            gene.fill_annotations(
                start=10 * (i - 1) + 1,
                stop=10 * i,
                strand="+",
                position=i,
                genetic_code=4,
            )
            genes.add(gene)
        yield genes

    def test_get_number_of_genes(self, family, genes):
        """Tests that the number of genes can be retrieved"""
        for gene in genes:
            family.add(gene)
        assert isinstance(len(family), int)
        assert len(family) == len(genes)

    @pytest.fixture
    def organisms(self, genes) -> Generator[Set[Organism], None, None]:
        """Create a set of organisms fill with genes to test edges"""
        organisms = set()
        genes = list(genes)
        nb_organisms = randint(2, 10)
        nb_genes_per_organisms = len(genes) // nb_organisms
        idx_org = 1
        contig_counter = 0
        while idx_org < nb_organisms:
            organism = Organism(f"organism_{idx_org}")
            contig = Contig(contig_counter, f"contig_{idx_org}")
            contig_counter
            organism.add(contig)
            idx_genes = 0
            while idx_genes < nb_genes_per_organisms:
                gene = genes[(idx_org - 1) * nb_genes_per_organisms + idx_genes]
                gene.fill_parents(organism, contig)
                contig[gene.start] = gene
                idx_genes += 1
            organisms.add(organism)
            idx_org += 1
        # last family fill with all the gene left
        organism = Organism(f"organism_{idx_org}")
        contig = Contig(contig_counter, f"contig_{idx_org}")
        organism.add(contig)
        idx_genes = (idx_org - 1) * nb_genes_per_organisms
        while idx_genes < len(genes):
            gene = genes[idx_genes]
            gene.fill_parents(organism, contig)
            contig[gene.start] = gene
            idx_genes += 1
        organisms.add(organism)
        yield organisms

    def test_get_org_dict(self, family, genes, organisms):
        """Tests that all organisms and genes are retrieved as expected"""
        for gene in genes:
            family.add(gene)
        org_dict = family.get_org_dict()
        assert isinstance(org_dict, dict)
        assert all(isinstance(org, Organism) for org in org_dict.keys())
        assert all(
            isinstance(gene, Gene)
            for gene_set in org_dict.values()
            for gene in gene_set
        )
        assert set(org_dict.keys()) == organisms
        assert (
            set([gene for gene_set in org_dict.values() for gene in gene_set]) == genes
        )

    def test_get_org_dict_with_no_organism_fill_to_genes(self, family, genes):
        """Tests that if genes are not fill with organism an AttributeError is returned"""
        for gene in genes:
            family.add(gene)
        with pytest.raises(AttributeError):
            _ = family.get_org_dict()

    def test_organisms(self, family, organisms, genes):
        """Tests that all organisms are retrieved as expected"""
        for gene in genes:
            family.add(gene)
        assert set(family.organisms) == organisms

    def test_number_of_organism(self, family, organisms, genes):
        """Tests that the expected number of organisms is found"""
        for gene in genes:
            family.add(gene)
        assert isinstance(family.number_of_organisms, int)
        assert family.number_of_organisms == len(organisms)

    def test_get_genes_per_org(self, family, organisms, genes):
        """Tests that for a giver organism, all the genes are retrieved as expected"""
        for gene in genes:
            family.add(gene)
        for organism in organisms:
            assert set(family.get_genes_per_org(organism)) == set(organism.genes)

    def test_get_genes_per_org_if_org_not_in_family(self, family):
        """Test that a KeyError is generated if an organism not belonging to the family is given"""
        with pytest.raises(KeyError):
            org = Organism("organism")
            _ = set(family.get_genes_per_org(org))

    @pytest.fixture
    def families(self, genes) -> Generator[Set[GeneFamily], None, None]:
        """Create a set of gene families fill with genes to test edges"""
        families = set()
        genes = list(genes)
        nb_families = randint(2, 10)
        nb_genes_per_family = len(genes) // nb_families
        idx_fam = 1
        while idx_fam < nb_families:
            family = GeneFamily(idx_fam, f"family_{idx_fam}")
            idx_genes = 0
            while idx_genes < nb_genes_per_family:
                gene = genes[(idx_fam - 1) * nb_genes_per_family + idx_genes]
                family.add(gene)
                gene.family = family
                idx_genes += 1
            families.add(family)
            idx_fam += 1
        # last family fill with all the gene left
        family = GeneFamily(idx_fam, f"family_{idx_fam}")
        idx_genes = (idx_fam - 1) * nb_genes_per_family
        while idx_genes < len(genes):
            gene = genes[idx_genes]
            family.add(gene)
            gene.family = family
            idx_genes += 1
        families.add(family)
        yield families

    @pytest.fixture
    def edges(self, families, genes, organisms) -> Generator[Set[Edge], None, None]:
        """Create a set of edges fill with genes and gene families to test edges"""
        edges = {}
        pair_genes = filter(
            lambda x: x[0] != x[1] and x[0].organism == x[1].organism,
            combinations_with_replacement(genes, 2),
        )
        for pair in pair_genes:
            key = frozenset([pair[0].family, pair[1].family])
            edge = edges.get(key)
            if edge is None:
                edge = Edge(pair[0], pair[1])
                edges[key] = edge
            else:
                edge.add_genes(pair[0], pair[1])
            pair[0].family.set_edge(pair[1].family, edge)
            pair[1].family.set_edge(pair[0].family, edge)
        yield set(edges.values())

    def test_get_neighbors_of_gene_family(self, families, edges):
        """Tests get all the expected neighbor of the family in the graph"""
        for family in families:
            assert all(
                isinstance(neighbor, GeneFamily) for neighbor in family.neighbors
            )
            expected_neighbors = set(
                [edge.source for edge in edges if edge.target == family]
            ).union(set([edge.target for edge in edges if edge.source == family]))
            assert set(family.neighbors) == expected_neighbors

    def test_get_number_of_neighbors(self, families, edges):
        """Tests that the expected number of neighbors is found"""
        for family in families:
            expected_neighbors = set(
                [edge.source for edge in edges if edge.target == family]
            ).union(set([edge.target for edge in edges if edge.source == family]))
            assert isinstance(family.number_of_neighbors, int)
            assert family.number_of_neighbors == len(expected_neighbors)

    #  Tests that the edges of a GeneFamily object can be retrieved
    def test_get_edges_of_gene_family(self, families, edges):
        """Tests that all the edges belonging to the family are retrieved"""
        for family in families:
            expected_edges = set(
                [
                    edge
                    for edge in edges
                    if edge.source == family or edge.target == family
                ]
            )
            assert all(isinstance(edge, Edge) for edge in family.edges)
            assert set(family.edges) == expected_edges

    def test_get_number_of_edges(self, families, edges):
        """Tests that the expected number of edges is found"""
        for family in families:
            expected_edges = set(
                [
                    edge
                    for edge in edges
                    if edge.source == family or edge.target == family
                ]
            )
            assert isinstance(family.number_of_edges, int)
            assert family.number_of_neighbors == len(expected_edges)

    def test_add_spot_to_gene_family(self, family):
        """Tests that a Spot object can be added to a GeneFamily object"""
        spot = Spot(1)
        family.add_spot(spot)
        assert spot in family.spots

    def test_add_non_spot_as_spot_in_family(self, family):
        """Tests that a non-spot object cannot be added to Gene Family"""
        with pytest.raises(TypeError):
            family.add_spot(323)

    def test_add_module_to_gene_family(self, family):
        """Tests that a Module object can be added to a GeneFamily object"""
        module = Module(1)
        family.set_module(module)
        assert module == family.module

    def test_add_non_module_as_module_in_family(self, family):
        """Tests that a non-module object cannot be added to Gene Family"""
        with pytest.raises(TypeError):
            family.set_module(323)

    # TODO test mk_bitarray
