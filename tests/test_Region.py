#! /usr/bin/env python3
import re

import pytest
from typing import Generator, Set
from random import randint

from ppanggolin.region import Region, Spot
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.genome import Gene, Contig, Organism


@pytest.fixture
def genes() -> Generator[Set[Gene], None, None]:
    """Create a set of genes to fill gene families
    """
    genes = set()
    for i in range(0, randint(11, 20)):
        gene = Gene(f"gene_{str(i)}")
        gene.fill_annotations(start=10 * i + 1, stop=10 * (i + 1), strand='+', position=i, genetic_code=4)
        genes.add(gene)
    yield genes


@pytest.fixture
def families(genes) -> Generator[Set[GeneFamily], None, None]:
    """Create a set of gene families fill with genes to test edges
    """
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
            family.add_gene(gene)
            gene.family = family
            idx_genes += 1
        families.add(family)
        idx_fam += 1
    # last family fill with all the gene left
    family = GeneFamily(idx_fam, f"family_{idx_fam}")
    idx_genes = (idx_fam - 1) * nb_genes_per_family
    while idx_genes < len(genes):
        gene = genes[idx_genes]
        family.add_gene(gene)
        gene.family = family
        idx_genes += 1
    families.add(family)
    yield families


class TestRegion:
    """Tests for region class
    """
    attr_val = {'score': 0, 'starter': None, 'stopper': None}

    @pytest.fixture
    def region(self) -> Generator[Region, None, None]:
        """Generate a region object to test class
        """
        yield Region("RGP")

    def test_cstr(self, region: Region):
        assert isinstance(region, Region)
        assert region.name == "RGP"
        assert isinstance(region._genes_getter, dict)
        for attr, value in self.attr_val.items():
            assert region.__getattribute__(attr) == value

    def test_add_gene(self, region):
        """Tests that genes can be aadded to a region
        """
        gene = Gene('gene')
        gene.fill_annotations(start=0, stop=10, strand='+', position=0)
        region[0] = gene
        assert len(region._genes_getter) == 1
        assert region._genes_getter[0] == gene
        assert region.starter == gene
        assert region.stopper == gene
        assert gene.RGP == region

    def test_add_gene_not_is_instance_gene(self, region):
        """Test that adding object with instance not Gene return a TypeError
        """
        with pytest.raises(TypeError):
            region[0] = 0

    def test_add_genes_at_position_already_taken(self, region):
        """Test that adding genes with same position return a ValueError
        """
        gene = Gene('gene')
        gene.fill_annotations(start=0, stop=10, strand='+', position=0)
        region[0] = gene
        with pytest.raises(ValueError):
            gene = Gene('gene')
            gene.fill_annotations(start=4, stop=12, strand='-', position=0)
            region[0] = gene

    def test_add_genes_from_different_contigs(self, region):
        """Test that adding genes from different contigs return an Exception
        """
        gene1, gene2 = Gene('gene_1'), Gene('gene_2')
        gene1.fill_annotations(start=0, stop=10, strand='+', position=0)
        gene2.fill_annotations(start=11, stop=20, strand='+', position=1)
        gene1.fill_parents(None, Contig('contig_1'))
        region[0] = gene1
        gene2.fill_parents(None, Contig('contig_2'))
        with pytest.raises(Exception):
            region[1] = gene2

    def test_add_genes_from_different_organisms(self, region):
        """Test that adding genes from different organisms return an Exception
        """
        gene1, gene2 = Gene('gene_1'), Gene('gene_2')
        gene1.fill_annotations(start=0, stop=10, strand='+', position=0)
        gene2.fill_annotations(start=11, stop=20, strand='+', position=1)
        gene1.fill_parents(Organism("org_1"))
        region[0] = gene1
        gene2.fill_parents(Organism("org_2"))
        with pytest.raises(Exception):
            region[1] = gene2

    def test_get_genes(self, region):
        """Tests that genes can be retrieved from the region
        """
        gene = Gene('gene')
        gene.fill_annotations(start=0, stop=10, strand='+', position=0)
        region[0] = gene
        assert region[0] == gene

    def test_get_genes_with_position_not_in_region(self, region):
        with pytest.raises(KeyError):
            _ = region[randint(0, 20)]

    def test_del_gene(self, region):
        """Tests that genes can be deleted from the region
        """
        gene = Gene('gene')
        gene.fill_annotations(start=0, stop=10, strand='+', position=0)
        region[0] = gene
        assert region[0] == gene
        del region[0]
        assert 0 not in region._genes_getter

    def test_get_length(self, region):
        """Tests that the length of the region can be retrieved
        """
        gene1, gene2 = Gene('gene_1'), Gene('gene_2')
        gene1.fill_annotations(start=0, stop=10, strand='+', position=0)
        gene2.fill_annotations(start=11, stop=20, strand='+', position=1)
        region[0] = gene1
        region[1] = gene2
        assert region.length == 20

    def test_get_organism(self, region):
        """Tests that the organism linked to the region can be retrieved
        """
        gene = Gene('gene')
        gene.fill_annotations(start=0, stop=10, strand='+', position=0)
        gene.fill_parents(Organism("org"))
        region[0] = gene
        assert region.organism.name == 'org'

    def test_get_contig(self, region):
        """Tests that the contig linked to the region can be retrieved
        """
        gene = Gene('gene')
        gene.fill_annotations(start=0, stop=10, strand='+', position=0)
        gene.fill_parents(None, Contig("contig"))
        region[0] = gene
        assert region.contig.name == 'contig'

    def test_is_whole_contig_true(self, region):
        """Tests that the property is_whole_contig return True if region is same length as contig
        """
        starter, stopper = Gene('starter'), Gene('stopper')
        starter.fill_annotations(start=0, stop=10, strand='+', position=0)
        stopper.fill_annotations(start=11, stop=20, strand='+', position=1)
        contig = Contig("contig")
        contig[starter.start], contig[stopper.start] = starter, stopper
        starter.fill_parents(None, contig), stopper.fill_parents(None, contig)
        region[starter.position], region[stopper.position] = starter, stopper
        assert region.is_whole_contig is True

    def test_is_whole_contig_false(self, region):
        """Tests that the property is_whole_contig return False if region is not same length as contig
        """
        before, starter, stopper, after = Gene('before'), Gene('starter'), Gene('stopper'), Gene('after')
        before.fill_annotations(start=0, stop=10, strand='+', position=0)
        starter.fill_annotations(start=11, stop=20, strand='+', position=1)
        stopper.fill_annotations(start=21, stop=30, strand='+', position=2)
        after.fill_annotations(start=31, stop=40, strand='+', position=3)
        contig = Contig("contig")
        contig[before.start], contig[after.start] = before, after
        contig[starter.start], contig[stopper.start] = starter, stopper
        before.fill_parents(None, contig), after.fill_parents(None, contig)
        starter.fill_parents(None, contig), stopper.fill_parents(None, contig)
        region[starter.position], region[stopper.position] = starter, stopper
        assert region.is_whole_contig is False

    def test_is_contig_border_true(self, region):
        """Test that property is_contig_border return true if the region is bordering the contig
        """
        before, starter, stopper, after = Gene('before'), Gene('starter'), Gene('stopper'), Gene('after')
        before.fill_annotations(start=0, stop=10, strand='+', position=0)
        starter.fill_annotations(start=11, stop=20, strand='+', position=1)
        stopper.fill_annotations(start=21, stop=30, strand='+', position=2)
        after.fill_annotations(start=31, stop=40, strand='+', position=3)
        contig = Contig("contig")
        before.fill_parents(None, contig), after.fill_parents(None, contig)
        starter.fill_parents(None, contig), stopper.fill_parents(None, contig)
        # Test bordering right
        contig[before.start], contig[starter.start], contig[stopper.start] = before, starter, stopper
        region[starter.position], region[stopper.position] = starter, stopper
        assert region.is_contig_border is True
        # Test bordering left
        del contig._genes_position[before.position]
        del contig._genes_getter[before.start]
        contig[after.start] = after
        assert region.is_contig_border is True

    def test_is_contig_border_false(self, region):
        """Tests that the property is_contig_border return False if region is not bordering the contig
        """
        before, starter, stopper, after = Gene('before'), Gene('starter'), Gene('stopper'), Gene('after')
        before.fill_annotations(start=0, stop=10, strand='+', position=0)
        starter.fill_annotations(start=11, stop=20, strand='+', position=1)
        stopper.fill_annotations(start=21, stop=30, strand='+', position=2)
        after.fill_annotations(start=31, stop=40, strand='+', position=3)
        contig = Contig("contig")
        contig[before.start], contig[after.start] = before, after
        contig[starter.start], contig[stopper.start] = starter, stopper
        before.fill_parents(None, contig), after.fill_parents(None, contig)
        starter.fill_parents(None, contig), stopper.fill_parents(None, contig)
        region[starter.position], region[stopper.position] = starter, stopper
        assert region.is_contig_border is False

    def test_len(self, region, genes):
        for gene in genes:
            region[gene.position] = gene
        assert isinstance(len(region), int)
        assert len(region) == len(genes)

    def test_equality(self, genes):
        """Test equality between two regions
        """
        region_1, region_2, region_3 = Region("RGP_1"), Region("RGP_2"), Region("RGP_3")
        max_pos = max(genes, key=lambda gene: gene.position).position
        for gene in genes:
            region_1[gene.position] = gene
            region_2[gene.position] = gene
            region_3[max_pos - gene.position + 1] = gene
        assert region_1 == region_2
        assert region_1 == region_3

    def test_not_equal(self, region, genes):
        """Test difference between two regions
        """
        for gene in genes:
            region[gene.position] = gene
        assert region != Region("other_RGP")

    def test_equality_with_not_instance_region(self, region):
        """Test comparison between a region and another object raise a TypeError
        """
        with pytest.raises(TypeError):
            assert region == 4

    def test_get_gene_families(self, region, genes, families):
        """Tests that gene families can be retrieved from the region
        """
        for gene in genes:
            region[gene.position] = gene
        assert all(isinstance(family, GeneFamily) for family in region.families)
        assert set(region.families) == families

    def test_get_number_of_gene_families(self, region, genes, families):
        """Tests that gene families can be retrieved from the region
        """
        for gene in genes:
            region[gene.position] = gene
        assert isinstance(region.number_of_families(), int)
        assert region.number_of_families() == len(families)

    # def test_get_bordering_genes(self, region, genes):
    #     # TODO test multigenic
    #     contig = Contig("contig")
    #     for gene in genes:
    #         contig[gene.start] = gene
    #         gene.fill_parents(None, contig)
    #         region[gene.position] = gene
    #     min_gene, max_gene = min(genes, key=lambda gene: gene.position), max(genes, key=lambda gene: gene.position)
    #     assert region.get_bordering_genes(1, {}) == [[min_gene], [max_gene]]


class TestSpot:
    @pytest.fixture
    def spot(self) -> Generator[Spot, None, None]:
        """Generate a spot for test
        """
        yield Spot(0)

    def test_cstr(self, spot):
        assert spot.ID == 0
        assert isinstance(spot._region_getter, dict) and len(spot._region_getter) == 0
        assert isinstance(spot._uniqOrderedSet, dict) and len(spot._uniqOrderedSet) == 0
        assert isinstance(spot._uniqContent, dict) and len(spot._uniqContent) == 0

    def test_cstr_assert_error(self):
        with pytest.raises(TypeError):
            Spot("spot_0")

    def test_repr(self, spot):
        """Test that the canonical string representing a spot does not change
        """
        assert repr(spot) == "Spot 0 - #RGP: 0"

    def test_str(self, spot):
        """Test that the writing spot method does not change
        """
        assert str(spot) == "spot_0"

    @pytest.fixture
    def region(self) -> Generator[Region, None, None]:
        """Create a region for test
        """
        yield Region("RGP_0")

    def test_add_region(self, spot, region):
        """Tests that adding a Region object to the Spot object works as expected
        """
        spot[region.name] = region
        assert region == spot._region_getter[region.name]

    def test_add_not_instance_region(self, spot):
        with pytest.raises(TypeError):
            spot["spot"] = "spot"

    def test_add_different_region_with_same_name(self, spot):
        """Test that adding a new Region same name than another in the spot return a KeyError
        """
        region_1, region_2 = Region("RGP"), Region("RGP")
        gene_1, gene_2 = Gene("gene_1"), Gene("gene_2")
        gene_1.fill_annotations(start=0, stop=10, strand='+', position=0)
        gene_2.fill_annotations(start=0, stop=10, strand='+', position=0)
        gene_1.family, gene_2.family = GeneFamily(0, "Fam_0"), GeneFamily(1, "Fam_1")
        region_1[0], region_2[0] = gene_1, gene_2
        spot[region_1.name] = region_1
        with pytest.raises(KeyError):
            spot[region_2.name] = region_2

    def test_add_two_time_the_same_region(self, spot, region):
        """Test that adding a new Region same name than another in the spot return a KeyError
        """
        gene = Gene("gene")
        gene.fill_annotations(start=0, stop=10, strand='+', position=0)
        gene.family = GeneFamily(0, "Fam")
        region[0] = gene
        spot[region.name] = region
        assert region in spot._region_getter.values()
        spot[region.name] = region
        assert region in spot._region_getter.values()

    def test_get_region(self, spot, region):
        """Tests that getting the region in the Spot object works as expected
        """
        spot[region.name] = region
        assert spot[region.name] == region

    def test_get_region_not_in_spot(self, spot):
        with pytest.raises(KeyError):
            _ = spot["rgp"]

    def test_delete_region_in_spot(self, spot, region):
        spot[region.name] = region
        del spot[region.name]
        assert region.name not in spot._region_getter

    def test_len(self, spot, region):
        assert isinstance(len(spot), int)
        assert len(spot) == 0
        spot[region.name] = region
        assert len(spot) == 1

    @pytest.fixture
    def regions(self, genes):
        regions = set()
        genes = sorted(list(genes), key=lambda x: x.position)
        nb_regions = randint(2, len(genes))
        nb_genes_per_region = len(genes) // nb_regions
        idx_region = 1
        while idx_region < nb_regions:
            region = Region(f"RGP_{idx_region}")
            genes_counter = 0
            while genes_counter < nb_genes_per_region:
                gene = genes.pop(0)
                region[gene.position] = gene
                gene.RGP = region
                genes_counter += 1
            regions.add(region)
            idx_region += 1
        # last region fill with all the gene left
        region = Region(f"RGP_{idx_region}")
        while len(genes) > 0:
            gene = genes.pop(0)
            region[gene.position] = gene
            gene.RGP = region
        regions.add(region)
        yield regions

    def test_get_all_regions(self, spot, regions):
        """Tests that getting all the region in the spot works as expected
        """
        for region in regions:
            spot[region.name] = region
        assert len(spot) == len(regions)
        assert all(type(region) == Region for region in spot.regions)
        assert regions == set(spot.regions)

    def test_get_families(self, spot, regions, families):
        """Tests that getting the gene families in the Spot object works as expected
        """
        for region in regions:
            spot[region.name] = region
        assert set(spot.families) == families

    def test_number_of_families(self, spot, regions, families):
        """Tests that getting the number of families in the spot works as expected
        """
        for region in regions:
            spot[region.name] = region
        assert spot.number_of_families() == len(families)

    def test_add_spot_to_families(self, spot, regions, families):
        """Tests that adding spot to families works as expected
        """
        for region in regions:
            spot[region.name] = region
        spot.spot_2_families()
        assert all(set(family.spots) == {spot} for family in spot.families)

    @pytest.fixture
    def srgps(self, regions):
        """Create a random number of same rgp for all regions
        """
        srgps = set()
        for region in regions:
            nb_sim_rgp = randint(1, 3)
            for idx_sim_rgp in range(1, nb_sim_rgp + 1):
                sim_rgp = Region(f"s{region.name}.{idx_sim_rgp}")
                for gene in region.genes:
                    sim_rgp[gene.position] = gene
                srgps.add(sim_rgp)
        yield srgps

    def test_get_uniq_rgp_set(self, spot, regions, families, srgps):
        """Tests that getting identical rgp in the Spot object works as expected
        """
        for region in list(regions) + list(srgps):  # With lists provide sRGP to be key RGP in dict
            spot[region.name] = region
        assert len(spot) == len(regions) + len(srgps)
        uniq2rgp = spot.get_uniq_to_rgp()
        for region, sim_rgps in uniq2rgp.items():
            assert region in regions
            assert set(region.families) == set.union(*[set(srgp.families) for srgp in sim_rgps])

    def test_get_uniq_ordered_set(self, spot, regions, families, srgps):
        """Tests that getting the unique synteny in the Spot object works as expected
        """
        for region in list(regions) + list(srgps):  # With lists provide sRGP to be key RGP in dict
            spot[region.name] = region
        assert len(spot) == len(regions) + len(srgps)
        assert spot.get_uniq_ordered_set().issubset(regions)

    def test_get_uniq_content(self, spot, regions, families, srgps):
        """Tests that getting the unique RGP in the Spot object works as expected
        """
        for region in list(regions) + list(srgps):  # With lists provide sRGP to be key RGP in dict
            spot[region.name] = region
        assert len(spot) == len(regions) + len(srgps)
        assert spot.get_uniq_ordered_set().issubset(regions)