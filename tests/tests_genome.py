#! /usr/bin/env python3

import pytest
from random import choices, randint, sample
from typing import Generator, Set, Tuple, Union
from pathlib import Path

from ppanggolin.genome import Feature, Gene, RNA, Contig, Organism
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region


class TestFeature:
    """Tests Feature class
    """

    def test_creation(self):
        """Tests that 'Feature' object is created successfully with the given identifier
        """
        feature = Feature('test_id')
        assert feature.ID == 'test_id'
        assert not feature.is_fragment
        assert feature.type == ''
        assert feature.start is None
        assert feature.stop is None
        assert feature.strand is None
        assert feature.product is None
        assert feature.name is None
        assert feature.local_identifier is None
        assert feature.organism is None
        assert feature.contig is None
        assert feature.dna is None

    def test_create_feature_assertion_error(self):
        """Tests that a Feature object cannot be created with a non-string type identifier"""
        with pytest.raises(AssertionError):
            Feature(4)

    def test_create_feature_empty_identifier(self):
        """Tests that a Feature object cannot be created with an empty identifier"""
        with pytest.raises(ValueError):
            Feature('')

    def test_fill_annotations(self):
        """Tests that 'fill_annotations' method fills the attributes correctly
        """
        feature = Feature('test_id')
        feature.fill_annotations(1, 10, '+', 'gene_type', 'name', 'product', 'local_id')
        assert feature.start == 1
        assert feature.stop == 10
        assert feature.type == 'gene_type'
        assert feature.strand == '+'
        assert feature.product == 'product'
        assert feature.name == 'name'
        assert feature.local_identifier == 'local_id'

    def test_set_organism_valid_type(self):
        """Tests that organism setter sets organism with valid type
        """
        feature = Feature('test')
        organism = Organism('organism')
        feature.organism = organism
        assert feature.organism == organism

    def test_set_organism_invalid_type(self):
        """Tests that organism setter return TypeError if sets organism with invalid type
        """
        feature = Feature('test')
        with pytest.raises(TypeError):
            feature.organism = 4

    def test_set_contig_valid_type(self):
        """Tests that contig setter sets contig with valid type
        """
        feature = Feature('test')
        contig = Contig('contig')
        feature.contig = contig
        assert feature.contig == contig

    def test_set_contig_invalid_type(self):
        """Tests that contig setter return TypeError if sets contig with invalid type
        """
        feature = Feature('test')
        with pytest.raises(TypeError):
            feature.contig = 4

    def test_fill_parents(self):
        """Tests that 'fill_parents' method associates the object with the given organism and contig
        """
        organism = Organism('org_id')
        contig = Contig('contig_name')
        feature = Feature('test_id')
        feature.fill_parents(organism, contig)
        assert feature.organism == organism
        assert feature.contig == contig

    def test_add_dna(self):
        """Tests that 'add_dna' method adds the DNA sequence to the object successfully
        """
        feature = Feature('test_id')
        feature.add_dna('ATCG')
        assert feature.dna == 'ATCG'

    def test_fill_annotations_type_error(self):
        """Tests that 'fill_annotations' method raises a TypeError if attribute value is not with correct type
        """
        feature = Feature('test_id')
        with pytest.raises(TypeError):
            feature.fill_annotations('1', 10, '+', 'gene_type', 'name', 'product', 'local_id')
        with pytest.raises(TypeError):
            feature.fill_annotations(1, "10", '+', 'gene_type', 'name', 'product', 'local_id')
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, 4, 'gene_type', 'name', 'product', 'local_id')
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, "+", 4, 'name', 'product', 'local_id')
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, '+', 'gene_type', 4, 'product', 'local_id')
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, '+', 'gene_type', 'name', 4, 'local_id')
        with pytest.raises(TypeError):
            feature.fill_annotations(1, 10, '+', 'gene_type', 'name', 'product', 4)

    def test_fill_annotations_value_error(self):
        """Tests that 'fill_annotations' method raises a TypeError if strand is not '+' or '-'
        """
        feature = Feature('test_id')
        with pytest.raises(ValueError):
            feature.fill_annotations(1, 10, '4', 'gene_type', 'name', 'product', 'local_id')

    def test_add_dna_type_error(self):
        """Tests that 'add_dna' method raises a TypeError if the DNA sequence is not a string
        """
        feature = Feature('test_id')
        with pytest.raises(AssertionError):
            feature.add_dna(123)

    def test_length_start_or_stop_are_not_known(self):
        """Tests that length property raises ValueError when start is not known
        """
        with pytest.raises(ValueError):
            feature = Feature('test')
            feature.stop = 10
            _ = feature.length
        with pytest.raises(ValueError):
            feature = Feature('test')
            feature.start = 1
            _ = feature.length


class TestGene:
    """Tests Gene class
    """
    def test_create_gene_object(self):
        """Tests that a Gene object can be created with a valid gene_id
        """
        gene = Gene('gene1')
        assert gene.ID == 'gene1'

    def test_fill_annotations(self):
        """Tests that Gene annotations can be filled with valid parameters
        """
        gene = Gene('gene1')
        gene.fill_annotations(start=1, stop=10, strand='+', position=10, genetic_code=4)
        assert gene.position == 10
        assert gene.genetic_code == 4

    #  Tests that Gene annotations cannot be filled with invalid parameters
    def test_fill_annotations_invalid_parameters(self):
        gene = Gene('gene1')
        with pytest.raises(TypeError):
            gene.fill_annotations(start=1, stop=10, strand='+', position='10', genetic_code=4)
        with pytest.raises(TypeError):
            gene.fill_annotations(start=1, stop=10, strand='+', position=10, genetic_code="4")

    def test_add_protein(self):
        """Tests that a protein sequence can be added to a Gene object
        """
        gene = Gene('gene1')
        gene.add_protein('MVKLAVLALALAVLALALALAVLALALAVLALALAVLALALAVLALALAVLALALAVLALALAVLALALA')
        assert gene.protein == 'MVKLAVLALALAVLALALALAVLALALAVLALALAVLALALAVLALALAVLALALAVLALALAVLALALA'

    def test_add_protein_non_string(self):
        """Tests that a non-string protein sequence cannot be added to a Gene object
        """
        gene = Gene('gene1')
        with pytest.raises(TypeError):
            gene.add_protein(123)

    def test_set_family_valid_type(self):
        """Tests that family setter sets family with valid type
        """
        gene = Gene('gene1')
        family = GeneFamily(0, 'family')
        gene.family = family
        assert gene.family == family

    def test_set_family_invalid_type(self):
        """Tests that family setter return TypeError if sets family with invalid type
        """
        gene = Gene('gene1')
        with pytest.raises(TypeError):
            gene.family = 4

    def test_set_rgp_valid_type(self):
        """Tests that RGP setter sets family with valid type
        """
        gene = Gene('gene1')
        region = Region(0)
        gene.RGP = region
        assert gene.RGP == region

    def test_set_rgp_invalid_type(self):
        """Tests that family setter return TypeError if sets family with invalid type
        """
        gene = Gene('gene1')
        with pytest.raises(TypeError):
            gene.RGP = 4