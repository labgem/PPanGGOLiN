#! /usr/bin/env python3

import pytest
from typing import Generator, Tuple
import gmpy2

from ppanggolin.genome import Feature, Gene, RNA, Contig, Organism
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region


class TestFeature:
    """Tests Feature class"""

    @pytest.fixture
    def feature(self) -> Generator[Feature, None, None]:
        """Generate a basic feature for tests"""
        yield Feature("test_id")

    def test_creation(self, feature):
        """Tests that 'Feature' is created successfully with the given identifier"""
        assert feature.ID == "test_id"
        assert not feature.is_fragment
        assert feature.type == ""
        assert feature.start is None
        assert feature.stop is None
        assert feature.strand is None
        assert feature.product is None
        assert feature.name is None
        assert feature.local_identifier is None
        assert feature.organism is None
        assert feature.contig is None
        assert feature.dna is None

    def test_create_feature_with_identifier_not_instance_string(self):
        """Tests that a Feature object cannot be created with a non-string type identifier"""
        with pytest.raises(AssertionError):
            Feature(4)

    def test_create_feature_empty_identifier(self):
        """Tests that a Feature object cannot be created with an empty identifier"""
        with pytest.raises(ValueError):
            Feature("")

    def tests_write_organism(self, feature):
        """Tests that write feature return feature name as string"""
        assert str(feature) == "test_id"

    def test_fill_annotations(self, feature):
        """Tests that 'fill_annotations' method fills the attributes correctly"""
        feature.fill_annotations(1, 10, "+", "gene_type", "name", "product", "local_id")
        assert feature.start == 1
        assert feature.stop == 10
        assert feature.type == "gene_type"
        assert feature.strand == "+"
        assert feature.product == "product"
        assert feature.name == "name"
        assert feature.local_identifier == "local_id"

    def test_fill_annotations_type_error(self, feature):
        """Tests that 'fill_annotations' method raises a TypeError if attribute value is not with the correct type"""
        with pytest.raises(TypeError):
            feature.fill_annotations(
                "1", 10, "+", "gene_type", "name", "product", "local_id"
            )
        with pytest.raises(TypeError):
            feature.fill_annotations(
                1, "10", "+", "gene_type", "name", "product", "local_id"
            )
        with pytest.raises(TypeError):
            feature.fill_annotations(
                1, 10, 4, "gene_type", "name", "product", "local_id"
            )
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, "+", 4, "name", "product", "local_id")
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, "+", "gene_type", 4, "product", "local_id")
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, "+", "gene_type", "name", 4, "local_id")
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, "+", "gene_type", "name", "product", 4)

    def test_fill_annotations_value_error(self, feature):
        """Tests that 'fill_annotations' method raises a TypeError if strand is not '+' or '-'"""
        with pytest.raises(ValueError):
            feature.fill_annotations(
                1, 10, "4", "gene_type", "name", "product", "local_id"
            )

    def test_fill_parents(self, feature):
        """Tests that 'fill_parents' method associates the object with the given organism and contig"""
        organism = Organism("org_id")
        contig = Contig(0, "contig_name")
        feature.fill_annotations(1, 10, "+", "gene_type", "name", "product", "local_id")
        feature.fill_parents(organism, contig)
        assert feature.organism == organism
        assert feature.contig == contig

    def test_fill_parents_with_organism_or_contig_only(self, feature):
        """Tests that Gene can be filled with only an organism or a contig"""
        organism = Organism("org")
        contig = Contig(0, "ctg")
        feature.fill_annotations(1, 10, "+", "gene_type", "name", "product", "local_id")
        feature.fill_parents(organism=organism)
        assert feature.organism == organism
        feature.fill_parents(contig=contig)
        assert feature.contig == contig

    def test_fill_parents_with_nothing(self, feature):
        """Tests that Gene cannot be filled with neither an organism and a contig"""
        with pytest.raises(AssertionError):
            feature.fill_parents()

    def test_set_organism(self, feature):
        """Tests that organism setter sets organism with the valid type"""
        organism = Organism("organism")
        feature.organism = organism
        assert feature.organism == organism

    def test_set_organism_not_isinstance_organism(self, feature):
        """Tests that organism setter return TypeError if sets organism with the invalid type"""
        with pytest.raises(TypeError):
            feature.organism = 4

    def test_set_contig(self, feature):
        """Tests that contig setter sets contig with the valid type"""
        contig = Contig(0, "contig")
        feature.contig = contig
        assert feature.contig == contig

    def test_set_contig_not_isinstance_contig(self, feature):
        """Tests that contig setter return TypeError if sets contig with the invalid type"""
        with pytest.raises(TypeError):
            feature.contig = 4

    def test_add_dna(self, feature):
        """Tests that 'add_dna' method adds the DNA sequence to the object successfully"""
        feature.add_sequence("ATCG")
        assert feature.dna == "ATCG"

    def test_add_dna_type_error(self, feature):
        """Tests that 'add_dna' method raises a TypeError if the DNA sequence is not a string"""
        with pytest.raises(AssertionError):
            feature.add_sequence(123)

    def test_length(self, feature):
        """Tests len method"""
        feature.fill_annotations(1, 10, "+", "gene_type", "name", "product", "local_id")
        assert isinstance(len(feature), int)
        assert len(feature) == 10

    def test_length_start_or_stop_are_not_known(self):
        """Tests that len raises ValueError when start is not known"""
        with pytest.raises(ValueError):
            feature = Feature("test")
            feature.stop = 10
            len(feature)
        with pytest.raises(ValueError):
            feature = Feature("test")
            feature.start = 1
            len(feature)

    @pytest.mark.parametrize(
        "coordinates, expected_overlaps_contig_edge_flag",
        [
            ([(1, 4), (3, 10)], False),
            ([(2, 4), (1, 1)], True),
            ([(1, 4), (1, 10)], False),
            ([(1, 4), (6, 10), (1, 2)], True),
            ([(5, 10), (9, 10), (1, 4)], True),
        ],
    )
    def test_overlaps_contig_edge(
        self, coordinates, expected_overlaps_contig_edge_flag
    ):
        feature = Feature("ID")
        feature.fill_annotations(start=1, stop=10, strand="+", coordinates=coordinates)

        assert feature.overlaps_contig_edge == expected_overlaps_contig_edge_flag


class TestRNA:
    """Tests RNA Class"""

    @pytest.fixture
    def rna(self) -> Generator[RNA, None, None]:
        """Generate a basic gene for tests"""
        yield RNA("rna")

    def test_create_gene_object(self, rna):
        """Tests that a Gene object can be created with a valid gene_id"""
        assert rna.ID == "rna"


class TestGene:
    """Tests Gene class"""

    @pytest.fixture
    def gene(self) -> Generator[Gene, None, None]:
        """Generate a basic gene for tests"""
        yield Gene("gene")

    def test_create_gene_object(self, gene):
        """Tests that a Gene object can be created with a valid gene_id"""
        assert gene.ID == "gene"
        assert gene.position is None
        assert gene._family is None
        assert gene._RGP is None
        assert gene.genetic_code is None
        assert gene.protein is None
        assert gene.is_partial is False
        assert gene._frame is None

    def test_fill_annotations(self, gene):
        """Tests that Gene annotations can be filled with valid parameters"""
        gene.fill_annotations(start=1, stop=10, strand="+", position=10, genetic_code=4)
        assert gene.position == 10
        assert gene.genetic_code == 4

    def test_fill_annotations_type_error(self, gene):
        """Tests that Gene annotations cannot be filled with invalid parameters"""
        with pytest.raises(TypeError):
            gene.fill_annotations(
                start=1, stop=10, strand="+", position="10", genetic_code=4
            )
        with pytest.raises(TypeError):
            gene.fill_annotations(
                start=1, stop=10, strand="+", position=10, genetic_code="4"
            )

    @pytest.mark.parametrize("frame", [0, 1, 2])
    def test_set_frame(self, frame):
        """Tests that frame can be set"""
        gene = Gene("gene")
        gene.frame = frame
        assert gene._frame == frame

    @pytest.mark.parametrize("frame", [0, 1, 2])
    def test_get_frame(self, frame):
        """Tests that frame can be getting"""
        gene = Gene("gene")
        gene.frame = frame
        assert gene.frame == frame

    def test_raise_assertion_error_if_frame_not_set(self):
        """Tests that frame cannot be return if it has not been set"""
        gene = Gene("gene")
        with pytest.raises(AssertionError):
            _ = gene.frame

    def test_raise_assertion_error_if_frame_already_set(self):
        """Tests that frame cannot be set if it has already been set"""
        gene = Gene("gene")
        gene.frame = 1
        with pytest.raises(AssertionError):
            gene.frame = 2

    @pytest.mark.parametrize("frame", [3, "1", 1.5])
    def test_raise_value_error_if_frame_not_0_1_or_2(self, frame):
        """Tests that frame cannot be set with value different from 0, 1 or 2"""
        gene = Gene("gene")
        with pytest.raises(ValueError):
            gene.frame = frame

    @pytest.mark.parametrize("frame", [0, 1, 2])
    def test_fill_partial_gene(self, frame):
        """Tests that Gene annotations can be filled with partial genes"""
        gene = Gene("gene")
        gene.fill_annotations(
            start=1, stop=10, strand="+", is_partial=True, frame=frame
        )
        assert gene.is_partial is True
        assert gene.frame == frame

    def test_add_protein(self, gene):
        """Tests that a protein sequence can be added to a Gene object"""
        gene.add_protein(
            "MVKLAVLALALAVLALALALAVLALALAVLALALAVLALALAVLALALAVLALALAVLALALAVLALALA"
        )
        assert (
            gene.protein
            == "MVKLAVLALALAVLALALALAVLALALAVLALALAVLALALAVLALALAVLALALAVLALALAVLALALA"
        )

    def test_add_protein_non_string(self, gene):
        """Tests that a non-string protein sequence cannot be added to a Gene object"""
        with pytest.raises(TypeError):
            gene.add_protein(123)

    def test_set_family(self, gene):
        """Tests that family setter sets family with the valid type"""
        family = GeneFamily(0, "family")
        gene.family = family
        assert gene.family == family

    def test_set_family_not_instance_gene_family(self, gene):
        """Tests that family setter return TypeError if sets family is not instance GeneFamily"""
        with pytest.raises(TypeError):
            gene.family = 4

    def test_set_rgp(self, gene):
        """Tests that RGP setter sets family with the valid type"""
        region = Region(0)
        gene.RGP = region
        assert gene.RGP == region

    def test_set_rgp_not_instance_region(self, gene):
        """Tests that family setter return TypeError if sets rgp is not instance Region"""
        with pytest.raises(TypeError):
            gene.RGP = 4


class TestContig:
    """Tests Contig class"""

    @pytest.fixture
    def contig(self) -> Generator[Contig, None, None]:
        """Generate basic contig for tests"""
        yield Contig(0, "contig")

    @pytest.fixture
    def gene(self) -> Generator[Gene, None, None]:
        """Generate basic gene for tests"""
        gene = Gene("test_gene")
        gene.fill_annotations(start=1, stop=10, strand="+", position=0, genetic_code=4)
        yield gene

    @pytest.fixture
    def genes(self) -> Generator[Tuple[Gene, Gene, Gene], None, None]:
        """Generate three basic genes for tests"""
        gene1 = Gene("test_gene1")
        gene1.fill_annotations(start=1, stop=10, strand="+", position=0, genetic_code=4)
        gene2 = Gene("test_gene2")
        gene2.fill_annotations(
            start=11, stop=20, strand="+", position=1, genetic_code=4
        )
        gene3 = Gene("test_gene3")
        gene3.fill_annotations(
            start=21, stop=30, strand="+", position=2, genetic_code=4
        )
        yield gene1, gene2, gene3

    def test_create_contig(self, contig):
        """Tests that a contig is correctly created"""
        assert contig.name == "contig"
        assert not contig.is_circular
        assert (
            contig._rna_getter == set()
        )  # Saving the rna annotations. We're not using them in the vast majority of cases.
        assert contig._genes_getter == {}
        assert contig._genes_position == []
        assert contig._organism is None

    def tests_write_contig(self, contig):
        """Tests that write contig return contig name as string"""
        assert str(contig) == "contig"

    def test_add_gene(self, gene, contig):
        """Tests that a gene can be added to the contig"""
        contig.add(gene)
        assert len(contig._genes_getter) == 1
        assert len(contig._genes_position) == 1
        assert contig._genes_getter[(gene.start, gene.stop, gene.strand)] == gene
        assert contig._genes_position[0] == gene

    def test_add_gene_at_far_position(self, gene, contig):
        """Tests that a gene can be added at each position and between position are fill with None"""
        contig.add(gene)
        new_gene = Gene("Gene2")
        new_gene.fill_annotations(
            start=50, stop=72, strand="+", position=6, genetic_code=4
        )
        contig.add(new_gene)
        assert len(contig._genes_position) == 7
        assert contig._genes_position[1:6] == [None] * 5

    def test_add_gene_not_instance_gene(self, contig):
        """Tests that the contig cannot be fill with a non-gene object"""
        with pytest.raises(TypeError):
            contig.add(1)
        with pytest.raises(TypeError):
            contig[1] = "4"

    def test_add_gene_with_start_already_taken(self, contig):
        """Tests that the contig cannot be fill with a non-gene object"""
        initial_gene = Gene("test_gene")
        initial_gene.fill_annotations(
            start=1, stop=12, strand="+", position=4, genetic_code=4
        )
        contig.add(initial_gene)
        with pytest.raises(ValueError):
            new_identical_gene = Gene("test_gene")
            new_identical_gene.fill_annotations(
                start=1, stop=12, strand="+", position=2, genetic_code=4
            )
            contig.add(new_identical_gene)

    def test_add_gene_without_position(self, contig):
        """Test that adding a gene not fill with position raise an AttributeError"""
        with pytest.raises(AttributeError):
            gene = Gene("test_gene")
            contig.add(gene)

    def test_number_of_genes(self, genes, contig):
        """Tests len method"""
        gene1, gene2, gene3 = genes
        contig.add(gene1)
        contig.add(gene2)
        contig.add(gene3)
        assert isinstance(contig.number_of_genes, int)
        assert contig.number_of_genes == 3

    def test_get_gene(self, gene, contig):
        """Tests that a gene can be retrieved by its position"""
        contig.add(gene)
        assert contig[0] == gene

    def test_get_genes(self, genes, contig):
        """Tests that a list of genes within a range can be retrieved"""
        gene1, gene2, gene3 = genes
        contig.add(gene1)
        contig.add(gene2)
        contig.add(gene3)
        assert set(contig.get_genes(0, 2)) == set(genes)

    def test_get_gene_with_non_integer_index(self, contig):
        """Tests that a gene cannot be retrieved with an index that is not an integer"""
        with pytest.raises(TypeError):
            _ = contig["a"]

    def test_get_genes_with_non_integer_begin_and_end_positions(self, genes, contig):
        """Tests that genes cannot be retrieved with non-integer begin and end positions"""
        gene1, gene2, gene3 = genes
        contig.add(gene1)
        contig.add(gene2)
        contig.add(gene3)
        with pytest.raises(TypeError):
            contig.get_genes("a", 2)
        with pytest.raises(TypeError):
            contig.get_genes(5, "b")
        with pytest.raises(TypeError):
            contig.get_genes("a", "b")

    def test_get_genes_with_end_position_lower_than_begin_position(self, genes, contig):
        """Tests that genes cannot be retrieved with end position lower than begin position"""
        gene1, gene2, gene3 = genes
        contig.add(gene1)
        contig.add(gene2)
        contig.add(gene3)
        with pytest.raises(ValueError):
            contig.get_genes(2, 0)

    def test_get_genes_with_end_position_greater_than_last_position(
        self, genes, contig
    ):
        """Tests that genes cannot be retrieved with given end position greater than last gene position in the contig"""
        gene1, gene2, gene3 = genes
        contig.add(gene1)
        contig.add(gene2)
        contig.add(gene3)
        with pytest.raises(IndexError):
            contig.get_genes(0, 3)

    def test_get_genes_with_end_position_greater_than_last_position_with_outrange_ok(
        self, genes, contig
    ):
        gene1, gene2, gene3 = genes
        contig.add(gene1)
        contig.add(gene2)
        contig.add(gene3)
        assert set(contig.get_genes(0, 5, outrange_ok=True)) == set(genes)

    def test_iterate_over_genes(self, genes, contig):
        """Tests that all genes in the contig can be iterated over"""
        gene1, gene2, gene3 = genes
        contig.add(gene1)
        contig.add(gene2)
        contig.add(gene3)
        assert list(contig.genes) == sorted(
            [gene1, gene2, gene3], key=lambda x: x.position
        )

    def test_add_rna(self, contig):
        """Tests that an RNA can be added to the contig"""
        rna = RNA("test_rna")
        contig.add_rna(rna)
        assert list(contig.RNAs) == [rna]

    def test_set_organism(self, contig):
        """Tests that an organism can be set to the contig"""
        organism = Organism("organism")
        contig.organism = organism
        assert contig.organism == organism

    def test_set_organism_with_not_instance_organism(self, contig):
        """Tests that the contig cannot be fill with a non-organism object"""
        with pytest.raises(TypeError):
            contig.organism = 4


class TestOrganism:
    """Tests Contig class"""

    @pytest.fixture
    def organism(self) -> Generator[Organism, None, None]:
        """Generate a basic organism for test"""
        yield Organism("organism")

    @pytest.fixture
    def contig(self) -> Generator[Contig, None, None]:
        """Generate a basic contig for test"""
        yield Contig(0, "contig")

    @pytest.fixture
    def gene(self) -> Generator[Gene, None, None]:
        """Generate a basic gene for test"""
        gene = Gene("test_gene")
        gene.fill_annotations(start=1, stop=10, strand="+", position=0, genetic_code=4)
        yield gene

    def test_create_organism(self, organism):
        """Tests that an Organism instance can be created with a valid name"""
        assert organism.name == "organism"
        assert organism._contigs_getter == {}
        assert organism._families is None
        assert organism.bitarray is None

    def test_create_organism_empty_name(self):
        """Tests that an Organism instance cannot be created with an empty name"""
        with pytest.raises(AssertionError):
            Organism("")

    def test_create_organism_with_name_not_string(self):
        """Tests that an Organism instance cannot be created with a name not instance string"""
        with pytest.raises(AssertionError):
            Organism(4)

    def tests_write_organism(self, organism):
        """Tests that write organism return organism name as string"""
        assert str(organism) == "organism"

    def test_add_contig(self, organism, contig):
        """Tests that a contig can be added to an Organism instance"""
        organism.add(contig)
        assert organism._contigs_getter["contig"] == contig

    def test_add_contig_not_instance_contig(self, organism):
        """Tests that a non Contig object cannot be added to an Organism instance"""
        with pytest.raises(AssertionError):
            organism.add(4)

    def test_add_contig_existing_name(self, organism, contig):
        """Tests that a contig with an existing name cannot be added to an Organism instance"""
        organism.add(contig)
        with pytest.raises(KeyError):
            organism.add(Contig(0, "contig"))

    def test_get_contig(self, organism, contig):
        """Tests that a contig can be retrieved from an Organism instance"""
        organism.add(contig)
        assert organism.get("contig") == contig

    def test_get_contig_not_instance_string(self, organism):
        """Tests that a non Contig object cannot be added to an Organism instance"""
        with pytest.raises(TypeError):
            organism.get(4)

    def test_get_nonexistent_contig(self, organism):
        """Tests that a non-existent contig cannot be retrieved from an Organism instance"""
        with pytest.raises(KeyError):
            organism.get("contig1")

    def test_number_of_contigs(self, organism):
        """Tests that the number of contigs in an organism instance can be retrieved"""
        organism.add(Contig(1, "contig1"))
        organism.add(Contig(2, "contig2"))

        assert organism.number_of_contigs == 2
        assert isinstance(len(organism), int)
        assert len(organism) == 2

    def test_get_families(self, organism, contig, gene):
        """Tests that gene families in an organism can be retrieved"""
        family = GeneFamily(0, "fam")
        family.add(gene)
        gene.fill_parents(organism, contig)
        organism.add(contig)
        contig[gene.start] = gene
        assert set(organism.families) == {family}

    def test_number_of_families(self, organism, contig, gene):
        """Tests that the number of gene families in an organism instance can be retrieved"""
        family = GeneFamily(0, "fam")
        family.add(gene)
        gene.fill_parents(organism, contig)
        organism.add(contig)
        contig.add(gene)
        assert organism.number_of_families() == 1

    def tests_get_genes(self, organism, contig, gene):
        """Tests that genes in an organism can be retrieved"""
        gene.fill_parents(organism, contig)
        organism.add(contig)
        contig.add(gene)
        assert set(organism.genes) == {gene}

    def test_number_of_genes(self, organism, contig, gene):
        """Tests that the number of genes in an organism instance can be retrieved"""
        gene.fill_parents(organism, contig)
        organism.add(contig)
        contig.add(gene)
        assert organism.number_of_genes() == 1

    def test_mk_bitarray(self, organism, contig):
        """Tests that a bitarray can be created for an Organism instance"""
        fam1 = GeneFamily(1, "fam1")
        fam2 = GeneFamily(2, "fam2")
        gene1 = Gene("gene1")
        gene2 = Gene("gene2")
        gene1.fill_annotations(start=1, stop=10, strand="+", position=0, genetic_code=4)
        gene2.fill_annotations(
            start=11, stop=19, strand="+", position=1, genetic_code=4
        )
        fam1.add(gene1)
        fam2.add(gene2)
        contig[gene1.start] = gene1
        contig[gene2.start] = gene2
        organism.add(contig)
        index = {fam1: 1, fam2: 2}
        organism.mk_bitarray(index)
        assert organism.bitarray == gmpy2.xmpz(6)
