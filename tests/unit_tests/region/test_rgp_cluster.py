#! /usr/bin/env python3

import pytest
from pyroaring import BitMap

from ppanggolin.RGP.rgp_cluster import RGPClustering, RegionProxy, RGPInfo, Contig



@pytest.fixture
def rgp_a():
    return RegionProxy(
        ID=1,
        name="RGP_A",
        families=BitMap([1, 2, 3, 4, 5]),
        is_contig_border=True,
        is_whole_contig=False,
    )


@pytest.fixture
def rgp_b():
    return RegionProxy(
        ID=2,
        name="RGP_B",
        families=BitMap([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        is_contig_border=False,
        is_whole_contig=False,
    )

@pytest.fixture
def clustering(tmp_path):
    clustering = RGPClustering(tmp_path / "pangenome.h5")
    clustering._fam_to_modules = {}
    clustering._rgp_to_nb_genes = {}
    clustering._rgp_contig_to_info = {
        "contig_1": Contig(
            organism="organism_1",
            is_circular=False,
            idx=0,
            name="contig_1",
        )
    }
    return clustering

@pytest.mark.parametrize(
    ("families_a", "families_b", "operation", "expected"),
    [
        (
            BitMap([1, 2, 3, 4, 5]),
            BitMap([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
            min,
            1.0,
        ),
        (
            BitMap([1, 2, 3, 4, 5]),
            BitMap([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
            max,
            0.5,
        ),
    ],
)
def test_grr(clustering, families_a, families_b, operation, expected):
    assert clustering._grr(families_a, families_b, operation) == expected



def test_rgp_metric(clustering, rgp_a, rgp_b):
    metric = clustering._rgp_metric(
        rgp_a,
        rgp_b,
        grr_cutoff=0,
        metric="min_grr",
    )

    assert metric is not None
    assert metric.min_grr == 1.0
    assert metric.max_grr == 0.5
    assert metric.incomplete_aware_grr == 1.0
    assert metric.shared_family == 5

def test_rgp_metric_below_cutoff(clustering, rgp_a, rgp_b):
    metric = clustering._rgp_metric(
        rgp_a,
        rgp_b,
        grr_cutoff=1.0,
        metric="max_grr",
    )

    assert metric is None


def test_rgp_metric_for_complete_rgps(clustering, rgp_a, rgp_b):
    rgp_a.is_contig_border = False
    rgp_a.is_whole_contig = False

    rgp_b.is_contig_border = False
    rgp_b.is_whole_contig = False

    metric = clustering._rgp_metric(
        rgp_a,
        rgp_b,
        grr_cutoff=0,
        metric="max_grr",
    )

    assert metric is not None
    assert metric.min_grr == 1.0
    assert metric.max_grr == 0.5
    assert metric.incomplete_aware_grr == 0.5
    assert metric.shared_family == 5


def test_construct_single(clustering):
    rgp = RGPInfo(
        name="RGP_1",
        families_ids=BitMap([1, 2, 3]),
        families=set(),
        contig="contig_1",
        is_contig_border=True,
        is_whole_contig=False,
    )

    clustering._rgp_to_nb_genes["RGP_1"] = 10

    region = clustering._construct_single(1, rgp)

    assert region.ID == 1
    assert region.name == "RGP_1"
    assert region.families == BitMap([1, 2, 3])
    assert region.organism == "organism_1"
    assert region.contig == "contig_1"
    assert region.length == 10
    assert region.is_contig_border
    assert not region.is_whole_contig
    assert not region.is_identical_region