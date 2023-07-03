#! /usr/bin/env python3

import pytest
from random import choices, randint, sample
from typing import Generator, Set
from pathlib import Path

from ppanggolin.genome import Gene, Organism, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.edge import Edge
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot, Module
from ppanggolin.metadata import Metadata


class TestPangenome:
    @pytest.fixture
    def pangenome(self) -> Generator[Pangenome, None, None]:
        """Create a pangenomes object for test

        :return: Generator with pangenomes object
        """
        pangenome = Pangenome()
        yield pangenome

    def test_cstr(self, pangenome):
        pangenome_attr_type = {
            "file": type(None),
            "_famGetter": dict,
            "_org_index": type(None),
            "_fam_index": type(None),
            "max_fam_id": int,
            "_orgGetter": dict,
            "_edgeGetter": dict,
            "_regionGetter": dict,
            "_spotGetter": dict,
            "_moduleGetter": dict,
            "status": dict,
            "parameters": dict
        }
        status_keys = [
            'genomesAnnotated',
            'geneSequences',
            'genesClustered',
            'defragmented',
            'geneFamilySequences',
            'neighborsGraph',
            'partitioned',
            'predictedRGP',
            'spots',
            'modules',
            "metadata",
            "metasources"
        ]
        metadata_keys = [
            "families",
            "genes",
            "genomes",
            "RGPs",
            "spots",
            "modules"
        ]
        for attr, attr_type in pangenome_attr_type.items():
            assert hasattr(pangenome, attr)
            assert isinstance(pangenome.__getattribute__(attr), attr_type)
            if attr_type == dict:
                if attr == "status":
                    assert len(pangenome.status) == len(status_keys)
                else:
                    assert len(pangenome.__getattribute__(attr)) == 0

        for status_key in status_keys:
            assert status_key in pangenome.status
            if status_key not in ["metadata", "metasources"]:
                assert pangenome.status[status_key] == "No"
            else:
                assert_res = "No" if status_key == "metadata" else []
                for metadata_key in metadata_keys:
                    assert metadata_key in pangenome.status[status_key]
                    assert pangenome.status[status_key][metadata_key] == assert_res
        assert pangenome.max_fam_id == 0

    def test_is_instance_pangenome(self, pangenome):
        assert isinstance(pangenome, Pangenome)


class TestPangenomeOrganism(TestPangenome):
    """Organism test"""
    def test_add_organism(self, pangenome):
        org = Organism("org")
        pangenome.add_organism(org)
        assert set(pangenome.organisms) == {org}

    def test_get_organism(self, pangenome):
        org = Organism("org")
        pangenome.add_organism(org)
        get_org = pangenome.get_organism("org")
        assert isinstance(get_org, Organism)
        assert org == get_org
        with pytest.raises(KeyError):
            pangenome.get_organism('organism')

    @pytest.fixture
    def orgs(self) -> Generator[Set[Organism], None, None]:
        """Create a list of organism object for test

        :return: Generator with list of organism object
        """
        orgs = set()
        for i in range(randint(5, 20)):
            org = Organism(str(i))
            orgs.add(org)
        yield orgs

    @pytest.fixture
    def add_organisms(self, pangenome, orgs):
        for org in orgs:
            pangenome.add_organism(org)

    def test_number_of_organisms(self, add_organisms, pangenome, orgs):
        assert isinstance(pangenome.number_of_organisms(), int)
        assert pangenome.number_of_organisms() == len(orgs)

    def test_add_organisms(self, add_organisms, pangenome, orgs):
        # 'set' because order is not guaranted
        # and org should be unique
        assert set(pangenome.organisms) == set(orgs)


class TestPangenomeGeneFamilies(TestPangenome):
    def test_max_fam_id_is_instance_int_and_egal_zero(self, pangenome):
        assert isinstance(pangenome.max_fam_id, int)
        assert pangenome.max_fam_id == 0

    def test_add_gene_family(self, pangenome):
        family = GeneFamily(pangenome.max_fam_id, "family")
        pangenome.add_gene_family(family)
        assert 1 == pangenome.max_fam_id
        with pytest.raises(KeyError):
            pangenome.add_gene_family(family)

    def test_get_gene_family(self, pangenome):
        family = GeneFamily(pangenome.max_fam_id, "family")
        pangenome.add_gene_family(family)
        assert isinstance(pangenome.get_gene_family("family"), GeneFamily)
        assert pangenome.get_gene_family("family") == family

    @pytest.fixture
    def families(self) -> Generator[Set[GeneFamily], None, None]:
        """Create a list of organism object for test

        :return: Generator with list of organism object
        """
        families = set()
        for i in range(randint(5, 20)):
            family = GeneFamily(family_id=i, name=f'family{i}')
            families.add(family)
        yield families

    @pytest.fixture
    def add_families(self, pangenome, families):
        for family in families:
            pangenome.add_gene_family(family)

    def test_number_of_gene_families_empty(self, add_families, pangenome, families):
        assert pangenome.number_of_gene_families() == len(families)


class TestPangenomeGene(TestPangenome):
    @pytest.fixture
    def genes(self):
        genes = set()
        for i in range(randint(5, 20)):
            gene = Gene(gene_id=i)
            genes.add(gene)
        yield genes

    def test_get_gene_empty(self, pangenome):
        with pytest.raises(KeyError):
            pangenome.get_gene(33)

    @pytest.fixture(name="organism_genes")
    def fill_org_with_genes(self):
        genes = set()
        organism = Organism(name="organism")
        for contig_id in range(randint(2, 10)):
            contig = organism.get_contig("k_{}".format(contig_id))
            for gene_idx in range(randint(2, 10)):
                gene = Gene(gene_id=f"{organism.name}.{contig_id}.{gene_idx}")
                gene.position = gene_idx
                gene.start = gene_idx
                contig.add_gene(gene)
                genes.add(gene)
        yield organism, genes

    @pytest.fixture(name="family_genes")
    def fill_family_with_genes(self, pangenome):
        genes = set()
        family = GeneFamily(family_id=pangenome.max_fam_id, name="family")
        for gene_idx in range(randint(2, 10)):
            gene = Gene(gene_id=f"{family.name}_{gene_idx}")
            gene.position = gene_idx
            gene.start = gene_idx
            family.add_gene(gene)
            genes.add(gene)
        yield family, genes

    def test_genes_organism_generator(self, pangenome, organism_genes):
        # orgs with genes.
        organism, genes = organism_genes
        pangenome.add_organism(organism)
        assert len(genes.difference(set(pangenome.genes))) == 0

    def test_get_gene_with_organism(self, pangenome, organism_genes):
        organism, genes = organism_genes
        pangenome.add_organism(organism)
        for gene in genes:
            assert pangenome.get_gene(gene.ID) == gene

    def test_genes_gene_families(self, family_genes, pangenome):
        """Genes are added in pan through their family."""
        family, genes = family_genes
        pangenome.add_gene_family(family)
        assert len(genes.difference(set(pangenome.genes))) == 0

    def test_get_with_gene_family(self, pangenome, family_genes):
        family, genes = family_genes
        pangenome.add_gene_family(family)
        for gene in genes:
            assert pangenome.get_gene(gene.ID) == gene

    def test_number_of_gene(self, pangenome, organism_genes):
        # orgs with genes.
        organism, genes = organism_genes
        pangenome.add_organism(organism)
        assert isinstance(pangenome.number_of_genes(), int)
        assert pangenome.number_of_genes() == len(genes)

    def test_get_multigenic(self, pangenome):
        # TODO make a better test
        multigenic = pangenome.get_multigenics(0.5)
        assert isinstance(multigenic, set)


class TestPangenomeEdge(TestPangenome):
    @staticmethod
    def make_gene_pair(gene_id_1: int = 1, gene_id_2: int = 2):
        """create a pair of genes that belong to the same organism in 2 different families."""
        gene1 = Gene(gene_id=f"gene_{gene_id_1}")
        gene2 = Gene(gene_id=f"gene_{gene_id_2}")
        fam1 = GeneFamily(family_id=1, name=f"fam_{gene_id_1}")
        fam2 = GeneFamily(family_id=2, name=f"fam_{gene_id_2}")
        ctg1 = Contig(name=f"ctg_{gene_id_1}")
        ctg2 = Contig(name=f"ctg_{gene_id_2}")
        fam1.add_gene(gene1)
        fam2.add_gene(gene2)
        organism = Organism(name=f"org_{choices([gene_id_1, gene_id_2], k=1)}")
        gene1.fill_parents(organism, ctg1)
        gene2.fill_parents(organism, ctg2)
        return gene1, gene2

    @pytest.fixture
    def gene_pair(self):
        return self.make_gene_pair()

    def test_add_edge(self, pangenome, gene_pair):
        gene1, gene2 = gene_pair
        edge = pangenome.add_edge(gene1, gene2)
        assert isinstance(edge, Edge)
        # addEdge doesn't act the same when the edge already exists.
        assert pangenome.add_edge(gene1, gene2) == edge

    def test_number_of_edges(self, pangenome, gene_pair):
        gene1, gene2 = gene_pair
        pangenome.add_edge(gene1, gene2)
        assert isinstance(pangenome.number_of_edges(), int)
        assert pangenome.number_of_edges() == 1

    def test_edges_one(self, pangenome, gene_pair):
        gene_1, gene_2 = gene_pair

        edges = []
        for _ in range(randint(2, 5)):
            edges.append(pangenome.add_edge(gene_1, gene_2))

        # always the same family couple
        #    = one edge, with several couple of genes
        # I use set because edges are uniques, it is not a multigraph.
        assert set(pangenome.edges) == set(edges)
        assert pangenome.number_of_edges() == 1

        edge = list(pangenome.edges).pop()
        assert edge.gene_pairs[0] == (gene_1, gene_2)

    @pytest.fixture
    def gene_pairs(self):
        gene_pairs = set()
        for _ in range(randint(5, 20)):
            gene_id_1, gene_id_2 = choices(range(randint(2, 10)), k=2)
            gene1, gene2 = self.make_gene_pair(gene_id_1, gene_id_2)
            gene_pairs.add((gene1, gene2))
        yield gene_pairs

    def test_edges_many_rand(self, pangenome, gene_pairs):
        edges = set()
        for gene_pair in gene_pairs:
            edges.add(pangenome.add_edge(*gene_pair))
        # I use set because edges are uniques, it is not a supergraph.
        assert set(pangenome.edges) == edges


class TestPangenomeBinary(TestPangenomeOrganism, TestPangenomeGeneFamilies):
    def test_get_org_index(self, add_organisms, pangenome, orgs):
        orgs_index = pangenome.get_org_index()
        assert isinstance(orgs_index, dict)
        index_know = set()
        for org, index in orgs_index.items():
            assert isinstance(org, Organism)
            assert isinstance(index, int)
            assert index not in index_know
            index_know.add(index)

    def test_compute_family_bitarrays_without_index_already_computed(self, add_families, pangenome):
        pangenome.compute_family_bitarrays()
        for family in pangenome.gene_families:
            assert family.bitarray is not None

    def test_compute_family_bitarrays_with_index_already_computed(self, add_families, pangenome):
        org_idx = pangenome.get_org_index()
        assert pangenome.compute_family_bitarrays() == org_idx


class TestPangenomeRGP(TestPangenome):
    def test_add_region(self, pangenome):
        rgp = Region(region_id="rgp")
        pangenome.add_region(rgp)
        assert len(pangenome._regionGetter) == 1
        assert pangenome._regionGetter["rgp"] == rgp

    def test_add_region_already_in_pangenome(self, pangenome):
        rgp = Region(region_id="rgp")
        pangenome.add_region(rgp)
        with pytest.raises(KeyError):
            pangenome.add_region(rgp)

    def test_get_region(self, pangenome):
        rgp = Region(region_id="rgp")
        pangenome.add_region(rgp)
        assert pangenome.get_region("rgp") == rgp

    def test_get_region_not_in_pangenome(self, pangenome):
        with pytest.raises(KeyError):
            pangenome.get_region("rgp")

    def test_number_of_rgp(self, pangenome):
        rgp = Region(region_id="rgp")
        pangenome.add_region(rgp)
        assert isinstance(pangenome.number_of_rgp(), int)
        assert pangenome.number_of_rgp() == 1


class TestPangenomeSpot(TestPangenome):
    def test_add_spot(self, pangenome):
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        assert len(pangenome._spotGetter) == 1
        assert pangenome._spotGetter[0] == spot

    def test_add_spot_already_in_pangenome(self, pangenome):
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        with pytest.raises(KeyError):
            pangenome.add_spot(spot)

    def test_get_spot_with_int(self, pangenome):
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        assert pangenome.get_spot(0) == spot

    def test_get_spot_with_str(self, pangenome):
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        assert pangenome.get_spot("spot_0") == spot

    def test_get_spot_not_in_pangenome(self, pangenome):
        with pytest.raises(KeyError):
            pangenome.get_spot(0)

    def test_number_of_spots(self, pangenome):
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        assert isinstance(pangenome.number_of_spots(), int)
        assert pangenome.number_of_spots() == 1


class TestPangenomeModule(TestPangenome):
    def test_add_module(self, pangenome):
        module = Module(module_id=0)
        pangenome.add_module(module)
        assert len(pangenome._moduleGetter) == 1
        assert pangenome._moduleGetter[0] == module

    def test_add_module_already_in_pangenome(self, pangenome):
        module = Module(module_id=0)
        pangenome.add_module(module)
        with pytest.raises(KeyError):
            pangenome.add_module(module)

    def test_get_module_with_int(self, pangenome):
        module = Module(module_id=0)
        pangenome.add_module(module)
        assert pangenome.get_module(0) == module

    def test_get_module_with_str(self, pangenome):
        module = Module(module_id=0)
        pangenome.add_module(module)
        assert pangenome.get_module("module_0") == module

    def test_get_module_not_in_pangenome(self, pangenome):
        with pytest.raises(KeyError):
            pangenome.get_module(0)

    def test_number_of_modules(self, pangenome):
        module = Module(module_id=0)
        pangenome.add_module(module)
        assert isinstance(pangenome.number_of_modules(), int)
        assert pangenome.number_of_modules() == 1

class TestPangenomeMetadata(TestPangenome):
    @pytest.fixture
    def add_element_to_pangenome(self, pangenome):
        metadata = Metadata(source="source", attribute="attr")
        family = GeneFamily(family_id=pangenome.max_fam_id, name="Fam")
        family.add_metadata(source=metadata.source, metadata=metadata)
        pangenome.add_gene_family(family)
        org = Organism("Org")
        org.add_metadata(source=metadata.source, metadata=metadata)
        ctg = org.get_contig("Ctg")
        gene = Gene("Gene")
        gene.position, gene.start = (0, 0)
        gene.add_metadata(source=metadata.source, metadata=metadata)
        ctg.add_gene(gene)
        pangenome.add_organism(org)
        rgp = Region("RGP")
        rgp.add_metadata(source=metadata.source, metadata=metadata)
        pangenome.add_region(rgp)
        spot = Spot(0)
        spot.add_metadata(source=metadata.source, metadata=metadata)
        pangenome.add_spot(spot)
        module = Module(0)
        module.add_metadata(source=metadata.source, metadata=metadata)
        pangenome.add_module(module)

    def test_metadata_sources(self, add_element_to_pangenome, pangenome):
        for metatype in ["families", "genomes", "genes", "RGPs", "spots", "modules"]:
            assert isinstance(pangenome.metadata_sources(metatype), set)
            assert pangenome.metadata_sources(metatype) == {'source'}

    def test_metadata(self, add_element_to_pangenome, pangenome):
        for metatype in ["families", "genomes", "genes", "RGPs", "spots", "modules"]:
            for metadata_gen in pangenome.metadata(metatype):
                for metadata in metadata_gen:
                    assert isinstance(metadata, Metadata)
                    assert metadata.source == 'source'

    def test_get_elem_by_metadata(self, add_element_to_pangenome, pangenome):
        for metatype, expected_type in {"families": GeneFamily, "genomes": Organism, "genes": Gene, "RGPs": Region,
                                        "spots": Spot, "modules": Module}.items():
            for elem in pangenome.get_elem_by_metadata(metatype, attribute="attr"):
                assert isinstance(elem, expected_type)
                for metadata in elem.metadata:
                    assert isinstance(metadata, Metadata)
                    assert metadata.source == 'source'

    def test_get_elem_by_sources(self, add_element_to_pangenome, pangenome):
        for metatype, expected_type in {"families": GeneFamily, "genomes": Organism, "genes": Gene, "RGPs": Region,
                                        "spots": Spot, "modules": Module}.items():
            for elem in pangenome.get_elem_by_sources(source='source', metatype=metatype):
                assert isinstance(elem, expected_type)
                for metadata in elem.metadata:
                    assert isinstance(metadata, Metadata)
                    assert metadata.source == 'source'
