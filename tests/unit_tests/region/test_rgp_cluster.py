#! /usr/bin/env python3

import pytest
from random import randint
from typing import Generator, Set
from ppanggolin.RGP import rgp_cluster
from ppanggolin.genome import Gene, Contig, Organism
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region
from ppanggolin.RGP.rgp_cluster import IdenticalRegions


@pytest.fixture
def genes() -> Generator[Set[Gene], None, None]:
    """Create a set of genes to fill gene families"""
    organism = Organism("organism")
    contig = Contig(0, "contig")
    genes = []
    for i in range(randint(11, 20)):
        gene = Gene(f"gene_{str(i)}")
        gene.fill_annotations(
            start=10 * i + 1, stop=10 * (i + 1), strand="+", position=i, genetic_code=4
        )
        gene.fill_parents(organism, contig)
        contig[gene.start] = gene
        genes.append(gene)
    return genes


@pytest.fixture
def families(genes) -> Generator[Set[GeneFamily], None, None]:
    """Create a set of gene families fill with genes to test edges"""
    families = set()
    genes = list(genes)
    nb_families = randint(9, 20)
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
def identical_rgps(genes, families) -> Generator[Set[Region], None, None]:
    """Create a set of identical rgps"""
    identical_rgps = set()
    for i in range(1, randint(6, 21)):
        rgp = Region(f"RGP_{i}")
        # the three rgp have the gene content.
        # in terms of family they are identical
        for gene in genes:
            rgp[gene.position] = gene
        identical_rgps.add(rgp)
    yield identical_rgps


class TestIdenticalRegions:
    def test_init_with_valid_inputs(self, identical_rgps, families):
        """Tests that the IdenticalRegions object is initialized correctly with valid inputs."""
        is_contig_border = True
        identical_regions = IdenticalRegions(
            "IdenticalRegions", identical_rgps, families, is_contig_border
        )

        assert identical_regions.name == "IdenticalRegions"
        assert identical_regions.rgps == identical_rgps
        assert identical_regions.families == families
        assert identical_regions.is_contig_border == is_contig_border

    @pytest.mark.parametrize("wrong_type", ["string", 1, 0.8, list(), dict()])
    def test_init_with_identical_rgps_not_isintance_set(self, wrong_type, families):
        """Tests that the IdenticalRegions object cannot be initialized with a not instance set for identical_rgps."""
        with pytest.raises(TypeError):
            IdenticalRegions("IdenticalRegions", wrong_type, families, True)

    def test_init_with_rgp_is_not_instance_region_in_identical_rgps(
        self, identical_rgps, families
    ):
        """Tests that the IdenticalRegions object raise TypeError if one element is not instance Region."""
        with pytest.raises(TypeError):
            IdenticalRegions(
                "IdenticalRegions", identical_rgps.union({1}), families, True
            )

    def test_init_with_empty_identical_rgps(self, families):
        """Tests that the IdenticalRegions object cannot be initialized with an empty set of identical regions."""
        with pytest.raises(ValueError):
            IdenticalRegions("IdenticalRegions", set(), families, True)

    @pytest.mark.parametrize("wrong_type", ["string", 1, 0.8, list(), dict()])
    def test_init_with_families_not_isintance_set(self, wrong_type, identical_rgps):
        """Tests that the IdenticalRegions object cannot be initialized with a not instance set."""
        with pytest.raises(TypeError):
            IdenticalRegions("IdenticalRegions", identical_rgps, wrong_type, True)

    def test_init_with_family_is_not_instance_genefamilies_in_families(
        self, identical_rgps, families
    ):
        """Tests that the IdenticalRegions object raise TypeError if one element is not instance Region."""
        with pytest.raises(TypeError):
            IdenticalRegions(
                "IdenticalRegions", identical_rgps, families.union({1}), True
            )

    def test_init_with_empty_families(self, identical_rgps):
        """Tests that the IdenticalRegions object cannot be initialized with an empty set of identical regions."""
        with pytest.raises(ValueError):
            IdenticalRegions("IdenticalRegions", identical_rgps, set(), True)

    def test_eq_with_equal_identical_regions(self):
        """Tests that the __eq__ method returns Trueè when comparing two IdenticalRegions objects that have the same families,
        identical regions, and contig border status.
        """
        rgp1 = Region("RGP1")
        rgp2 = Region("RGP2")
        family1 = GeneFamily(1, "Family1")
        family2 = GeneFamily(2, "Family2")
        identical_rgps1 = {rgp1, rgp2}
        identical_rgps2 = {rgp1, rgp2}
        families1 = {family1, family2}
        families2 = {family1, family2}
        is_contig_border = True

        identical_regions1 = IdenticalRegions(
            "IdenticalRegions", identical_rgps1, families1, is_contig_border
        )
        identical_regions2 = IdenticalRegions(
            "IdenticalRegions", identical_rgps2, families2, is_contig_border
        )

        assert identical_regions1 == identical_regions2

    def test_eq_with_non_identical_regions(self):
        """Tests that the __eq__ method returns False when comparing
        two IdenticalRegions objects that have different families.
        """
        rgp1 = Region("RGP1")
        rgp2 = Region("RGP2")
        family1 = GeneFamily(1, "Family1")
        family2 = GeneFamily(2, "Family2")
        identical_rgps1 = {rgp1, rgp2}
        identical_rgps2 = {rgp1, rgp2}
        families1 = {family1, family2}
        families2 = {family1}
        is_contig_border = True

        identical_regions1 = IdenticalRegions(
            "IdenticalRegions", identical_rgps1, families1, is_contig_border
        )
        identical_regions2 = IdenticalRegions(
            "IdenticalRegions", identical_rgps2, families2, is_contig_border
        )

        assert identical_regions1 != identical_regions2


def test_compute_grr():
    """Tests that compute_grr returns the correct value when there is a non-zero intersection between families"""
    set1 = {1, 2, 3, 4, 5}
    set2 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

    assert rgp_cluster.compute_grr(set1, set2, min) == 1.0
    assert rgp_cluster.compute_grr(set1, set2, max) == 0.5


def test_dereplicate_rgp(identical_rgps):
    list_identical_rgps = list(identical_rgps)
    rgp1 = list_identical_rgps[0]
    assert rgp_cluster.dereplicate_rgp({rgp1}) == [rgp1]

    identical_region_obj = rgp_cluster.IdenticalRegions(
        name="identical_rgps_0",
        identical_rgps=identical_rgps,
        families=set(list_identical_rgps[0].families),
        is_contig_border=True,
    )
    assert rgp_cluster.dereplicate_rgp(rgps=identical_rgps)[0] == identical_region_obj


def test_compute_rgp_metric(genes, families):
    RGP_a = Region("A")
    RGP_b = Region("B")
    list_genes = sorted(genes, key=lambda x: x.position)

    for g in list_genes[:8]:
        RGP_a[g.position] = g
    for g in list_genes[3:7]:
        RGP_b[g.position] = g

    assert RGP_a.is_contig_border
    assert not RGP_b.is_contig_border

    shared_families = len(set(RGP_a.families).intersection(set(RGP_b.families)))
    expected_grr = (
        RGP_a.ID,
        RGP_b.ID,
        {
            "incomplete_aware_grr": shared_families
            / min(len(set(RGP_a.families)), len(set(RGP_b.families))),
            "min_grr": shared_families
            / min(len(set(RGP_a.families)), len(set(RGP_b.families))),
            "max_grr": shared_families
            / max(len(set(RGP_a.families)), len(set(RGP_b.families))),
            "shared_family": shared_families,
        },
    )
    # min_grr
    min_result = rgp_cluster.compute_rgp_metric(RGP_a, RGP_b, 0, "min_grr")
    assert min_result == expected_grr

    # incomplete_aware_grr: same as min grr as rgp1 is incomplete
    incomplete_aware_result = rgp_cluster.compute_rgp_metric(
        RGP_a, RGP_b, 0, "incomplete_aware_grr"
    )

    assert incomplete_aware_result == expected_grr

    # max grr is below cutoff so None is returned
    assert rgp_cluster.compute_rgp_metric(RGP_a, RGP_b, 1000, "max_grr") is None
