#! /usr/bin/env python3

import pytest
from random import choices, randint
from typing import Generator, Set, Tuple, Union

from ppanggolin.genome import Gene, Organism, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.edge import Edge
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot, Module
from ppanggolin.metadata import Metadata


class TestPangenome:
    """This class tests methods in pangenome class associated to pangenome direclty.
    For pangenome components, there are subclasses to test each component.
    This class also generate a pangenome for all the test
    """

    @pytest.fixture
    def pangenome(self) -> Generator[Pangenome, None, None]:
        """Create a pangenomes object for test

        :return: Generator with the pangenome object
        """
        pangenome = Pangenome()
        yield pangenome

    def test_cstr(self, pangenome):
        """
        Tests the constructor method of the Pangenome class.
        It checks that all attributes are present and have the correct type and or value.

        :param pangenome: Test the function

        :return: A pangenome object
        """
        pangenome_attr_type = {
            "file": type(None),
            "_fam_getter": dict,
            "_org_index": type(None),
            "_fam_index": type(None),
            "_max_fam_id": int,
            "_org_getter": dict,
            "_edge_getter": dict,
            "_region_getter": dict,
            "_spot_getter": dict,
            "_module_getter": dict,
            "status": dict,
            "parameters": dict,
        }
        status_keys = [
            "genomesAnnotated",
            "geneSequences",
            "genesClustered",
            "defragmented",
            "geneFamilySequences",
            "neighborsGraph",
            "partitioned",
            "predictedRGP",
            "spots",
            "modules",
            "metadata",
            "metasources",
        ]
        metadata_keys = ["families", "genes", "genomes", "RGPs", "spots", "modules"]
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
        """Tests whether the pangenome object is an instance of the Pangenome class.
        This test is important because it ensures that the class name does not change and that we are working
        with a Pangenome object, and not some other type of object.

        :param pangenome: Object to test if is an instance of the pangenome class

        :raise AssertionError: If pangenome is not an instance of the pangenome class
        """
        assert isinstance(pangenome, Pangenome)

    def test_add_file_is_not_path(self, pangenome):
        """Tests that the add_file method raises an AssertionError if a file is not an instance of the Path class

        :param pangenome: Pangenome object to test method
        """
        with pytest.raises(AssertionError):
            pangenome.add_file("pangenome.h5")


class TestPangenomeOrganism(TestPangenome):
    """This class tests methods in pangenome class associated to organisms."""

    @pytest.fixture
    def organism(self) -> Generator[Organism, None, None]:
        """Create a basic organism"""
        yield Organism(name="organism")

    def test_add_organism(self, pangenome, organism):
        """Tests the add_organism method of the Pangenome class.

        :param pangenome: Pangenome object to test method
        :param organism: organism object to test method
        """
        pangenome.add_organism(organism)
        assert set(pangenome.organisms) == {organism}

    def test_add_organism_already_in_pangenome(self, pangenome, organism):
        """Tests that adding organism that already exist return a KeyError.

        :param pangenome: Pangenome object to test method
        :param organism: organism object to test method
        """
        pangenome.add_organism(organism)
        with pytest.raises(KeyError):
            pangenome.add_organism(organism)

    def test_add_organism_not_instance_organism(self, pangenome):
        """Ensure that it raises an AssertionError when a non-Organism object is passed as an argument.

        :param pangenome: Pangenome object to test method
        """
        with pytest.raises(AssertionError):
            pangenome.add_organism("org")

    def test_get_organism(self, pangenome, organism):
        """Tests the get_organism method of the Pangenome class.

        :param pangenome: Pangenome object to test method
        :param organism: organism object to test method
        """
        pangenome.add_organism(organism)
        get_org = pangenome.get_organism("organism")
        assert isinstance(get_org, Organism)
        assert organism == get_org

    def test_get_organism_not_in_pangenome(self, pangenome):
        """Ensure that it raises a KeyError when an Organism is not in the pangenome.

        :param pangenome: Pangenome object to test method
        """
        with pytest.raises(KeyError):
            pangenome.get_organism("org")

    def test_get_organism_with_name_not_instance_string(self, pangenome):
        """Ensure that it raises an AssertionError when a non-string name is passed as organism name.

        :param pangenome: Pangenome object to test method
        """
        with pytest.raises(AssertionError):
            pangenome.get_organism(33)

    @pytest.fixture
    def organisms(self) -> Generator[Set[Organism], None, None]:
        """Create a set of organism object for test

        :return: Generator with the set of organism object
        """
        orgs = set()
        for i in range(randint(5, 20)):
            org = Organism(str(i))
            orgs.add(org)
        yield orgs

    @pytest.fixture
    def add_organisms(self, pangenome, organisms):
        """Add the set of organims to pangenome

        :param pangenome: Pangenome object to test method
        :param organisms: Set of organisms to add to pangenome
        """
        for org in organisms:
            pangenome.add_organism(org)

    def test_number_of_organisms(self, add_organisms, pangenome, organisms):
        """Tests the number_of_organisms method of the pangenome class.

        :param add_organisms: Method to add organisms
        :param pangenome: Pangenome object to test method
        :param organisms: Set of organisms to add to pangenome
        """
        assert isinstance(pangenome.number_of_organisms, int)
        assert pangenome.number_of_organisms == len(organisms)


class TestPangenomeGeneFamilies(TestPangenome):
    """This class tests methods in pangenome class associated to gene families."""

    @pytest.fixture
    def family(self) -> Generator[GeneFamily, None, None]:
        """Create a Gene Family object

        :return: Generator with a Gene Family object
        """
        family = GeneFamily(0, "family")
        yield family

    def test_max_fam_id_is_instance_int_and_egal_zero(self, pangenome):
        """Tests that the max_fam_id attribute is corretly set

        :param pangenome: Pangenome object to test method
        """
        assert isinstance(pangenome.max_fam_id, int)
        assert pangenome.max_fam_id == 0

    def test_add_gene_family(self, pangenome, family):
        """Tests the add_gene_family method of the Pangenome class.

        :param pangenome: Pangenome object to test method
        :param family: gene family object to test method
        """
        pangenome.add_gene_family(family)
        assert 1 == pangenome.max_fam_id
        assert set(pangenome.gene_families) == {family}

    def test_add_gene_family_already_in_pangenome(self, pangenome, family):
        """Tests that adding gene family that already exist return a KeyError.

        :param pangenome: Pangenome object to test method
        :param family: gene family object to test method
        """
        pangenome.add_gene_family(family)
        with pytest.raises(KeyError):
            pangenome.add_gene_family(family)

    def test_get_gene_family(self, pangenome, family):
        """Tests that get_gene_family return a gene family object corresponding to the requested gene family

        :param pangenome: Pangenome object to test method
        :param family: gene family object to test method
        """
        pangenome.add_gene_family(family)
        assert isinstance(pangenome.get_gene_family("family"), GeneFamily)
        assert pangenome.get_gene_family("family") == family

    def test_get_gene_family_not_in_pangenome(self, pangenome, family):
        """Tests that return a KeyError if family does not exist in pangenome

        :param pangenome: Pangenome object to test method
        :param family: gene family object to test method
        """
        with pytest.raises(KeyError):
            pangenome.get_gene_family("fam")

    def test_get_gene_family_with_name_not_isinstance_string(self, pangenome):
        """Tests that return an AssertionError if family name used to get family is not string

        :param pangenome: Pangenome object to test method
        """
        with pytest.raises(AssertionError):
            pangenome.get_gene_family(3)

    @pytest.fixture
    def families(self) -> Generator[Set[GeneFamily], None, None]:
        """Create a set of Gene Family object for test

        :return: Generator with the set of organism object
        """
        families = set()
        for i in range(randint(5, 20)):
            family = GeneFamily(family_id=i, name=f"family{i}")
            families.add(family)
        yield families

    @pytest.fixture
    def add_families(self, pangenome, families):
        """Add the set of gene families to pangenome

        :param pangenome: pangenome object to test method
        :param families: set of gene families to add to pangenome
        """
        for family in families:
            pangenome.add_gene_family(family)

    def test_number_of_gene_families_empty(self, add_families, pangenome, families):
        """Tests the number_of_gene_families method of the pangenome class.

        :param add_families: Method to add gene families
        :param pangenome: pangenome object to test method
        :param families: set of families to add to pangenome
        """
        assert isinstance(pangenome.number_of_gene_families, int)
        assert pangenome.number_of_gene_families == len(families)


class TestPangenomeGene(TestPangenome):
    """This class tests methods in pangenome class associated to Gene."""

    @pytest.fixture
    def genes(self) -> Generator[Set[Gene], None, None]:
        """Create a set of Gene object for test

        :return: Generator with the set of organism object
        """
        genes = set()
        for i in range(randint(5, 20)):
            gene = Gene(gene_id=i)
            genes.add(gene)
        yield genes

    @pytest.fixture(name="organism_genes")
    def fill_org_with_genes(self) -> Generator[Union[Organism, Set[Gene]], None, None]:
        """Fill an organism with a random set of gene

        :return: Organism with genes
        """
        genes = set()
        organism = Organism(name="organism")
        for ctg_counter, contig_id in enumerate(range(randint(2, 10))):
            contig = Contig(ctg_counter, "k_{}".format(contig_id))
            organism.add(contig)
            for gene_idx in range(randint(2, 10)):
                gene = Gene(gene_id=f"{organism.name}.{contig_id}.{gene_idx}")
                gene.position = gene_idx
                gene.start = gene_idx
                contig[gene.start] = gene
                genes.add(gene)
        yield organism, genes

    @pytest.fixture(name="family_genes")
    def fill_family_with_genes(self, pangenome):
        """Fill a gene family with a random set of gene

        :return: Gene family with genes
        """
        genes = set()
        family = GeneFamily(family_id=pangenome.max_fam_id, name="family")
        for gene_idx in range(randint(2, 10)):
            gene = Gene(gene_id=f"{family.name}_{gene_idx}")
            gene.position = gene_idx
            gene.start = gene_idx
            family.add(gene)
            genes.add(gene)
        yield family, genes

    def test_genes_generator_from_organism(self, pangenome, organism_genes):
        """Tests genes generator from organism in the pangenome object

        :param pangenome: Pangenome object
        :param organism_genes: method to get an organism object filled with genes
        """
        organism, genes = organism_genes
        pangenome.add_organism(organism)
        assert genes == set(pangenome.genes)

    def test_get_gene_with_organism(self, pangenome, organism_genes):
        """Tests get genes from organism in pangenome object

        :param pangenome: Pangenome object
        :param organism_genes: Method to get an organism object filled with genes
        """
        organism, genes = organism_genes
        pangenome.add_organism(organism)
        for gene in genes:
            assert pangenome.get_gene(gene.ID) == gene

    def test_genes_generator_from_gene_families(self, family_genes, pangenome):
        """Tests genes generator from gene families in pangenome object

        :param pangenome: Pangenome object to test method
        :param family_genes: method to get a gene family object filled with genes
        """
        family, genes = family_genes
        pangenome.add_gene_family(family)
        assert genes == set(pangenome.genes)

    def test_get_with_gene_family(self, pangenome, family_genes):
        """Tests genes generator from gene families in pangenome object

        :param pangenome: Pangenome object to test method
        :param family_genes: method to get a gene family object filled with genes
        """
        family, genes = family_genes
        pangenome.add_gene_family(family)
        for gene in genes:
            assert pangenome.get_gene(gene.ID) == gene

    def test_get_gene_not_in_pangenome(self, pangenome):
        """Tests that return a KeyError if gene does not exist in pangenome

        :param pangenome: Pangenome object to test method
        """
        with pytest.raises(KeyError):
            pangenome.get_gene("12151405613024")

    def test_get_gene_with_id_not_string(self, pangenome):
        """Tests that return an AssertionError if gene identifier is not a string

        :param pangenome: Pangenome object to test method
        """
        with pytest.raises(AssertionError):
            pangenome.get_gene(gene_id=4)

    def test_number_of_genes(self, pangenome, organism_genes):
        """Tests get number of genes in pangenome object

        :param pangenome: pangenome object to test method
        :param organism_genes: method to get a organism object fill with genes
        """
        organism, genes = organism_genes
        pangenome.add_organism(organism)
        assert isinstance(pangenome.number_of_genes, int)
        assert pangenome.number_of_genes == len(genes)

    def test_get_multigenic(self, pangenome):
        # TODO make a better test
        """Tests get multigenic genes in pangenome object

        :param pangenome: pangenome object to test method
        """
        multigenic = pangenome.get_multigenics(0.5)
        assert isinstance(multigenic, set)


class TestPangenomeEdge(TestPangenome):
    """This class tests methods in pangenome class associated to Edge."""

    @staticmethod
    def make_gene_pair(gene_id_1: int = 1, gene_id_2: int = 2) -> Tuple[Gene, Gene]:
        """Create a pair of genes that belong to the same organism in two different families

        :return: Two genes linked to contigs, organism and gene families
        """
        gene1 = Gene(gene_id=f"gene_{gene_id_1}")
        gene2 = Gene(gene_id=f"gene_{gene_id_2}")
        fam1 = GeneFamily(family_id=1, name=f"fam_{gene_id_1}")
        fam2 = GeneFamily(family_id=2, name=f"fam_{gene_id_2}")
        ctg1 = Contig(1, name=f"ctg_{gene_id_1}")
        ctg2 = Contig(2, name=f"ctg_{gene_id_2}")
        fam1.add(gene1)
        fam2.add(gene2)
        organism = Organism(name=f"org_{choices([gene_id_1, gene_id_2], k=1)}")
        gene1.fill_parents(organism, ctg1)
        gene2.fill_parents(organism, ctg2)
        return gene1, gene2

    @pytest.fixture
    def gene_pair(self) -> Generator[Tuple[Gene, Gene], None, None]:
        """Call method to create a pair of genes that belong to the same organism in two different families

        :return: Two genes linked to contigs, organism and gene families
        """
        yield self.make_gene_pair()

    def test_add_edge(self, pangenome, gene_pair):
        """Tests the add_edge method of the Pangenome class.

        :param pangenome: Pangenome object to test method
        :param gene_pair: Pair of gene coding for the edge
        """
        gene1, gene2 = gene_pair
        edge = pangenome.add_edge(gene1, gene2)
        assert isinstance(edge, Edge)
        assert set(pangenome.edges) == {edge}

    def test_add_edge_already_in_pangenome(self, pangenome, gene_pair):
        """Tests that adding the same pair of genes as edge return the edge.

        :param pangenome: Pangenome object to test method
        :param gene_pair: Pair of gene coding for the edge
        """
        gene1, gene2 = gene_pair
        edge = pangenome.add_edge(gene1, gene2)
        assert pangenome.add_edge(gene1, gene2) == edge

    def test_add_edge_with_gene_not_isinstance_gene(self, pangenome):
        """Tests that return an AssertionError if genes are not Gene objects

        :param pangenome: Pangenome object to test method
        """
        with pytest.raises(AssertionError):
            pangenome.add_edge("gene1", "gene2")

    def test_number_of_edges(self, pangenome, gene_pair):
        """Tests the number_of_edges method of the Pangenome class.

        :param pangenome: Pangenome object to test method
        :param gene_pair: Pair of gene coding for the edge
        """
        pangenome.add_edge(*gene_pair)
        assert isinstance(pangenome.number_of_edges, int)
        assert pangenome.number_of_edges == 1


class TestPangenomeBinary(TestPangenomeOrganism, TestPangenomeGeneFamilies):
    """This class tests methods in pangenome class associated to binary methods."""

    # TODO Better test for this part
    def test_get_org_index(self, add_organisms, pangenome):
        """Tests the get_org_index function in pangenome class

        :param add_organisms: Add organisms to the pangenome
        :param pangenome: Pass the pangenome object
        """
        orgs_index = pangenome.get_org_index()
        assert isinstance(orgs_index, dict)
        index_know = set()
        for org, index in orgs_index.items():
            assert isinstance(org, Organism)
            assert isinstance(index, int)
            assert index not in index_know
            index_know.add(index)

    def test_compute_family_bitarrays_with_index_already_computed(
        self, add_organisms, add_families, pangenome
    ):
        """Tests the compute_family_bitarrays function in Pangenome class

        :param add_families: Add families to the pangenome object
        :param pangenome: Access the pangenome object
        """
        org_idx = pangenome.get_org_index()
        assert pangenome.compute_family_bitarrays() == org_idx

    def test_compute_family_bitarrays_without_index_already_computed(
        self, add_organisms, add_families, pangenome
    ):
        """Tests the compute_family_bitarrays function of the Pangenome class.

        :param add_families: Add families to the pangenome
        :param pangenome: Test the compute_family_bitarrays function
        """
        pangenome.compute_family_bitarrays()
        for family in pangenome.gene_families:
            assert family.bitarray is not None

    def test_get_fam_index(self, add_families, pangenome):
        """Tests the get_org_index function in pangenome class

        :param add_families: Add families to the pangenome
        :param pangenome: Pass the pangenome object
        """
        fams_index = pangenome.get_fam_index()
        assert isinstance(fams_index, dict)
        index_know = set()
        for fam, index in fams_index.items():
            assert isinstance(fam, GeneFamily)
            assert isinstance(index, int)
            assert index not in index_know
            index_know.add(index)

    def test_compute_org_bitarrays_with_index_already_computed(
        self, add_organisms, add_families, pangenome
    ):
        """Tests the compute_family_bitarrays function in Pangenome class

        :param add_families: Add families to the pangenome object
        :param pangenome: Access the pangenome object
        """
        fams_index = pangenome.get_fam_index()
        assert pangenome.compute_org_bitarrays() == fams_index

    def test_compute_org_bitarrays_without_index_already_computed(
        self, add_organisms, add_families, pangenome
    ):
        """Tests the compute_family_bitarrays function of the Pangenome class.

        :param add_families: Add families to the pangenome
        :param pangenome: Test the compute_family_bitarrays function
        """
        pangenome.compute_org_bitarrays()
        for organism in pangenome.organisms:
            assert organism.bitarray is not None


class TestPangenomeRGP(TestPangenome):
    """This class tests methods in pangenome class associated to Region"""

    def test_add_region(self, pangenome):
        """Tests the add_region method in the Pangenome class.

        :param pangenome: Access the pangenome object
        """
        rgp = Region(name="rgp")
        pangenome.add_region(rgp)
        assert len(pangenome._region_getter) == 1
        assert pangenome._region_getter["rgp"] == rgp

    def test_add_region_already_in_pangenome(self, pangenome):
        """Tests that adding region already in pangenome return a KeyError.

        :param pangenome: Access the pangenome object
        """
        rgp = Region(name="rgp")
        pangenome.add_region(rgp)
        with pytest.raises(KeyError):
            pangenome.add_region(rgp)

    def test_add_region_with_isinstance_not_region(self, pangenome):
        """Tests that adding an object with not Region type return an AssertionError.

        :param pangenome: Access the pangenome object
        """
        with pytest.raises(AssertionError):
            pangenome.add_region("rgp")

    def test_get_region(self, pangenome):
        """Tests the get_region method in the Pangenome class.

        :param pangenome: Access the pangenome object
        """
        rgp = Region(name="rgp")
        pangenome.add_region(rgp)
        assert pangenome.get_region("rgp") == rgp

    def test_get_region_not_in_pangenome(self, pangenome):
        """Tests get the region not in pangenome return a KeyError.

        :param pangenome: Access the pangenome object
        """
        with pytest.raises(KeyError):
            pangenome.get_region("rgp")

    def test_get_region_with_isinstance_not_string(self, pangenome):
        """Tests that getting a region with not string as identifier return an AssertionError.

        :param pangenome: Access the pangenome object
        """
        with pytest.raises(AssertionError):
            pangenome.get_region(15646)

    def test_number_of_rgp(self, pangenome):
        """Tests the number_of_rgp method in the Pangenome class.

        :param pangenome: Pass the pangenome object to the function
        """
        rgp = Region(name="rgp")
        pangenome.add_region(rgp)
        assert isinstance(pangenome.number_of_rgp, int)
        assert pangenome.number_of_rgp == 1


class TestPangenomeSpot(TestPangenome):
    """This class tests methods in pangenome class associated to Spot."""

    def test_add_spot(self, pangenome):
        """Tests the add_spot method in the Pangenome class.

        :param pangenome: Access the pangenome object
        """
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        assert len(pangenome._spot_getter) == 1
        assert pangenome._spot_getter[0] == spot

    def test_add_spot_already_in_pangenome(self, pangenome):
        """Tests that adding spot already in pangenome return a KeyError.

        :param pangenome: Access the pangenome object
        """
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        with pytest.raises(KeyError):
            pangenome.add_spot(spot)

    def test_add_spot_with_isinstance_not_spot(self, pangenome):
        """Tests that adding an object with not Spot type return an AssertionError.

        :param pangenome: Access the pangenome object
        """
        with pytest.raises(AssertionError):
            pangenome.add_spot(4564)

    def test_get_spot_with_int(self, pangenome):
        """Tests get_spot method with integer in pangenome class

        :param pangenome: Access the pangenome object
        """
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        assert pangenome.get_spot(0) == spot

    def test_get_spot_with_str(self, pangenome):
        """Tests get_spot method with string in pangenome class

        :param pangenome: Access the pangenome object
        """
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        assert pangenome.get_spot("spot_0") == spot

    def test_get_spot_not_in_pangenome(self, pangenome):
        """Tests that getting spot not in pangenome return a KeyError.

        :param pangenome: Access the pangenome object
        """
        with pytest.raises(KeyError):
            pangenome.get_spot(544654)

    def test_number_of_spots(self, pangenome):
        """Tests number_of_spots methods in Pangenome class

        :param pangenome: Access the pangenome object
        """
        spot = Spot(spot_id=0)
        pangenome.add_spot(spot)
        assert isinstance(pangenome.number_of_spots, int)
        assert pangenome.number_of_spots == 1


class TestPangenomeModule(TestPangenome):
    """This class tests methods in pangenome class associated to Modules."""

    def test_add_module(self, pangenome):
        """Tests the add_module method in the Pangenome class.

        :param pangenome: Access the pangenome object
        """
        module = Module(module_id=0)
        pangenome.add_module(module)
        assert len(pangenome._module_getter) == 1
        assert pangenome._module_getter[0] == module

    def test_add_module_already_in_pangenome(self, pangenome):
        """Tests that adding module already in pangenome return a KeyError.

        :param pangenome: Access the pangenome object
        """
        module = Module(module_id=0)
        pangenome.add_module(module)
        with pytest.raises(KeyError):
            pangenome.add_module(module)

    def test_add_module_with_isinstance_not_region(self, pangenome):
        """Tests that adding an object with not Module type return an AssertionError.

        :param pangenome: Access the pangenome object
        """
        with pytest.raises(AssertionError):
            pangenome.add_module("module")

    def test_get_module_with_int(self, pangenome):
        """Tests get_module method with integer in pangenome class

        :param pangenome: Access the pangenome object
        """
        module = Module(module_id=0)
        pangenome.add_module(module)
        assert pangenome.get_module(0) == module

    def test_get_module_with_str(self, pangenome):
        """Tests get_module method with string in pangenome class

        :param pangenome: Access the pangenome object
        """
        module = Module(module_id=0)
        pangenome.add_module(module)
        assert pangenome.get_module("module_0") == module

    def test_get_module_not_in_pangenome(self, pangenome):
        """Tests that getting module not in pangenome return a KeyError.

        :param pangenome: Access the pangenome object
        """
        with pytest.raises(KeyError):
            pangenome.get_module(0)

    def test_number_of_modules(self, pangenome):
        """Tests number_of_modules methods in Pangenome class

        :param pangenome: Access the pangenome object
        """
        module = Module(module_id=0)
        pangenome.add_module(module)
        assert isinstance(pangenome.number_of_modules, int)
        assert pangenome.number_of_modules == 1


class TestPangenomeMetadata(TestPangenome):
    """This class tests methods in pangenome class associated to Metadata."""

    @pytest.fixture
    def add_element_to_pangenome(self, pangenome):
        """Adds a metadata element to each element of pangenome

        :param pangenome: Access the pangenome object
        """
        metadata = Metadata(source="source", attribute="attr")
        family = GeneFamily(family_id=pangenome.max_fam_id, name="Fam")
        family.add_metadata(metadata=metadata)
        pangenome.add_gene_family(family)
        org = Organism("Org")
        org.add_metadata(metadata=metadata)
        ctg = Contig(0, "Ctg")
        org.add(ctg)
        gene = Gene("Gene")
        gene.fill_annotations(start=1, stop=100, position=0, strand="+")
        gene.add_metadata(metadata=metadata)
        ctg.add(gene)
        pangenome.add_organism(org)
        rgp = Region("RGP")
        rgp.add_metadata(metadata=metadata)
        pangenome.add_region(rgp)
        spot = Spot(0)
        spot.add_metadata(metadata=metadata)
        pangenome.add_spot(spot)
        module = Module(0)
        module.add_metadata(metadata=metadata)
        pangenome.add_module(module)

    def test_select_elem(self, add_element_to_pangenome, pangenome):
        """Tests the select_elem method of the Pangenome class.

        :param add_element_to_pangenome: Add elements to the pangenome
        :param pangenome: Access the pangenome object
        """
        assert all(
            isinstance(elem, GeneFamily)
            for elem in set(pangenome.select_elem("families"))
        )
        assert all(
            isinstance(elem, Organism) for elem in set(pangenome.select_elem("genomes"))
        )
        assert all(
            isinstance(elem, Gene) for elem in set(pangenome.select_elem("genes"))
        )
        assert all(
            isinstance(elem, Region) for elem in set(pangenome.select_elem("RGPs"))
        )
        assert all(
            isinstance(elem, Spot) for elem in set(pangenome.select_elem("spots"))
        )
        assert all(
            isinstance(elem, Module) for elem in set(pangenome.select_elem("modules"))
        )
        with pytest.raises(KeyError):
            pangenome.select_elem("error")

    def test_metadata_sources(self, add_element_to_pangenome, pangenome):
        """Tests the metadata_sources method of the Pangenome class.

        :param add_element_to_pangenome: Add elements to the pangenome
        :param pangenome: Access the pangenome object
        """
        for metatype in ["families", "genomes", "genes", "RGPs", "spots", "modules"]:
            assert isinstance(pangenome.metadata_sources(metatype), set)
            assert pangenome.metadata_sources(metatype) == {"source"}

    def test_metadata(self, add_element_to_pangenome, pangenome):
        """Tests the metadata generator of the Pangenome class.

        :param add_element_to_pangenome: Add elements to the pangenome
        :param pangenome: Access the pangenome object
        """
        for metatype in ["families", "genomes", "genes", "RGPs", "spots", "modules"]:
            for metadata_gen in pangenome.metadata(metatype):
                for metadata in metadata_gen:
                    assert isinstance(metadata, Metadata)
                    assert metadata.source == "source"

    def test_get_elem_by_metadata(self, add_element_to_pangenome, pangenome):
        """Tests the metadata generator filtered by metadata attribute of the Pangenome class.

        :param add_element_to_pangenome: Add elements to the pangenome
        :param pangenome: Access the pangenome object
        """
        for metatype, expected_type in {
            "families": GeneFamily,
            "genomes": Organism,
            "genes": Gene,
            "RGPs": Region,
            "spots": Spot,
            "modules": Module,
        }.items():
            for elem in pangenome.get_elem_by_metadata(metatype, attribute="attr"):
                assert isinstance(elem, expected_type)
                for metadata in elem.metadata:
                    assert isinstance(metadata, Metadata)
                    assert metadata.source == "source"

    def test_get_elem_by_source(self, add_element_to_pangenome, pangenome):
        """Tests the metadata generator filtered by source of the Pangenome class.

        :param add_element_to_pangenome: Add elements to the pangenome
        :param pangenome: Access the pangenome object
        """
        for metatype, expected_type in {
            "families": GeneFamily,
            "genomes": Organism,
            "genes": Gene,
            "RGPs": Region,
            "spots": Spot,
            "modules": Module,
        }.items():
            for elem in pangenome.get_elem_by_source(
                source="source", metatype=metatype
            ):
                assert isinstance(elem, expected_type)
                for metadata in elem.metadata:
                    assert isinstance(metadata, Metadata)
                    assert metadata.source == "source"
