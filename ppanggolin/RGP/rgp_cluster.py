#!/usr/bin/env python3

# default libraries
import logging
import argparse
import os
from itertools import combinations
from collections.abc import Callable
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Union, Any
from pathlib import Path

# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd

# local libraries
from ppanggolin import RGP
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region, Spot, Module
from ppanggolin.formats import check_pangenome_info
from ppanggolin.utils import restricted_float, mk_outdir, Timer
from ppanggolin.geneFamily import GeneFamily

import typing as tp
from dataclasses import dataclass, field
from pyroaring import BitMap
from ppanggolin.formats.h5reader import H5Reader, TableAttribute

class RegionProxy:

    def __init__(
        self,
        ID: int,
        name: str,
        families: BitMap,
        is_contig_border: bool,
        is_whole_contig: bool,
        children=None,
        modules=None,
        contig=None,
        organism=None,
        length=0,
    ):
        self.ID: int = ID
        self.name: str = name
        self.families: BitMap = families
        self.children: set[RegionProxy] = children
        self.modules: set[int] = modules

        if self.children:
            self.organism: str = next(iter(self.children)).organism
            self.contig: str = next(iter(self.children)).contig
        else:
            self.organism: str = organism
            self.contig: str = contig

        self.length: int = length
        self.nb_families: int = len(self.families)

        self.is_contig_border: bool = is_contig_border
        self.is_whole_contig: bool = is_whole_contig

    @property
    def is_identical_region(self) -> bool:
        return self.children is not None and len(self.children) > 0

    def __repr__(self):
        return f"RegionProxy2(ID={self.ID}, name='{self.name}')"

    def __str__(self):
        return self.name

    def __hash__(self) -> int:
        return id(self)

    def __eq__(self, rhs: "RegionProxy") -> bool:
        return all(self.families == rhs.families,
                   self.children == rhs.children,
                   self.is_contig_border == rhs.is_contig_border)

    def __lt__(self, obj):
        return self.ID < obj.ID

    def __gt__(self, obj):
        return self.ID > obj.ID

    def __le__(self, obj):
        return self.ID <= obj.ID

    def __ge__(self, obj):
        return self.ID >= obj.ID

@dataclass
class RGPTable:
    rgp: tp.Annotated[str, TableAttribute(name="RGP", transform=lambda x: x.decode('utf-8'))]
    gene: tp.Annotated[str, TableAttribute(name="gene", transform=lambda x: x.decode('utf-8'))]
    _table: str = "/RGP"

@dataclass
class GeneFamTable:
    gene: tp.Annotated[str, TableAttribute(name="gene", transform=lambda x: x.decode('utf-8'))]
    family: tp.Annotated[str, TableAttribute(name="geneFam", transform=lambda x: x.decode('utf-8'))]
    _table: str = "/geneFamilies"

@dataclass
class AnnotationsGeneTable:
    name: tp.Annotated[str, TableAttribute(name="ID", transform=lambda x: x.decode('utf-8'))]
    contig: tp.Annotated[int, TableAttribute(name="contig")]
    _table: str = "/annotations/genes"

@dataclass
class RGPSpotTable:
    rgp: tp.Annotated[str, TableAttribute(name="RGP", transform=lambda x: x.decode('utf-8'))]
    spot: tp.Annotated[int, TableAttribute(name="spot")]
    _table: str = "/spots"

@dataclass
class ModuleTable:
    fam: tp.Annotated[str, TableAttribute(name="geneFam", transform=lambda x: x.decode('utf-8'))]
    module: tp.Annotated[int, TableAttribute(name="module")]
    _table: str = "/modules"

@dataclass
class ContigTable:
    genome: tp.Annotated[str, TableAttribute(name="genome", transform=lambda x: x.decode('utf-8'))]
    contig: tp.Annotated[str, TableAttribute(name="name", transform=lambda x: x.decode('utf-8'))]
    is_circular: tp.Annotated[bool, TableAttribute(name="is_circular")]
    idx: tp.Annotated[int, TableAttribute(name="ID")]
    _table: str = "/annotations/contigs"

@dataclass
class GenesTable:
    name: tp.Annotated[str, TableAttribute(name="ID", transform=lambda x: x.decode('utf-8'))]
    genedata: tp.Annotated[int, TableAttribute(name="genedata_id")]
    contig: tp.Annotated[int, TableAttribute(name="contig")]
    _table: str = "/annotations/genes"

@dataclass
class GeneDataTable:
    # gene_type: tp.Annotated[str, TableAttribute(name="gene_type", transform=lambda x: x.decode('utf-8'))]
    idx: tp.Annotated[int, TableAttribute(name="genedata_id")]
    start: tp.Annotated[int, TableAttribute(name="start")]
    stop: tp.Annotated[int, TableAttribute(name="stop")]
    position: tp.Annotated[int, TableAttribute(name="position")]
    _table: str = "/annotations/genedata"

@dataclass
class RGPGeneProxy:
    start: int
    stop: int
    position: int

@dataclass
class RGPGenes:
    contig: int
    is_circular_contig: bool   
    genes: list[RGPGeneProxy]


@dataclass
class RGPInfo:
    name: str
    families: set[str]
    families_ids: BitMap
    is_contig_border: bool
    is_whole_contig: bool
    contig: str


@dataclass
class ContigBorderPosition:
    last_gene_position: int
    last_gene_idx: int
    first_gene_idx: int
    first_gene_position: int
    gene_count: int


@dataclass
class ContigBorderGenes:
    first_gene: str
    last_gene: str
    gene_count: int


@dataclass
class RGPMetric:
    max_grr: float
    min_grr: float
    incomplete_aware_grr: float
    shared_family: int


@dataclass
class Contig:
    organism: str
    is_circular: bool
    idx: int
    name: str


RGPMetricType = tp.Literal["max_grr", "min_grr", "incomplete_aware_grr"]

@dataclass
class RGPClusteringOptions:
    grr_cutoff: float = 0.3
    metric: RGPMetricType = "incomplete_aware_grr"
    unmerge_identical_rgps: bool = True
    output: Path = Path("rgp_clustering")
    basename: str = "rgp_cluster"
    graph_formats: list[str] = ("gexf", "graphml")
    with_metadata: bool = False
    metadata_sources: list[str] = field(default_factory=list)
    metadata_sep: str = "|"

class RGPClustering:
    def __init__(self, pangenome_h5: Path):
        self.h5 = Path(pangenome_h5)
        self.rgps: set[RegionProxy] = set()
        self.metrics: list[RGPMetric] = []
        self.identical_regions: int = 0
        self.graph: nx.Graph = nx.Graph()
        self._rgp_to_spot: dict[str, int] = None
        self._fam_to_modules: dict[str, set[str]] = None
        self._contig_to_organism: dict[str, str] = None
        self._rgp_to_nb_genes: dict[str, int] = None
        self._rgp_to_contig_info: dict[str, RGPGenes] = None

    def _get_rgp_spot(self, reader: H5Reader) -> dict[str, int]:
        rgp_to_spot: dict[str, int] = {}
        for table in reader.fetch(RGPSpotTable):
            rgp_to_spot[table.rgp] = table.spot
        return rgp_to_spot

    def _get_fam_to_modules(self, reader: H5Reader) -> dict[str, set[str]]:
        fam_to_modules: dict[str, set[str]] = defaultdict(set)
        for table in reader.fetch(ModuleTable):
            fam_to_modules[table.fam].add(table.module)
        return fam_to_modules

    def _get_contig_to_organism(self, reader: H5Reader) -> dict[str, str]:
        contig_to_organism: dict[str, str] = {}
        for table in reader.fetch(ContigTable):
            contig_to_organism[table.contig] = table.genome
        return contig_to_organism

    def _get_contig_to_is_circular(self, reader: H5Reader) -> dict[str, bool]:
        circular_contig_ids: dict[str, bool] = {}
        for table in reader.fetch(ContigTable):
            circular_contig_ids[table.idx] = table.is_circular
        return circular_contig_ids

    def _get_contig_to_info(
        self, reader: H5Reader, contigs_to_keep: set[str] = None
    ) -> dict[str, Contig]:

        contig_to_info: dict[str, Contig] = {}
        for table in reader.fetch(ContigTable):
            if contigs_to_keep is None or table.contig in contigs_to_keep:
                contig_to_info[table.contig] = Contig(
                    organism=table.genome,
                    is_circular=table.is_circular,
                    idx=table.idx,
                    name=table.contig,
                )

        return contig_to_info

    def _get_contig_border_genes(
        self, reader: H5Reader
    ) -> dict[str, ContigBorderPosition]:
        """
        Get contig information such as if it is circular, last gene position and last gene idx.
        """

        contig_ids_with_rgp = {
            contig_info.idx for contig_info in self._rgp_contig_to_info.values()
        }
        genedata_to_contig_ids: dict[int, list[int]] = defaultdict(list)

        # contig_genedata_id_to_gene_name: dict[tuple[int, int], str] = {}

        # Map gene name with genetadata id
        # and contig id with genmetadata id in annotations/genes table

        for table in reader.fetch(GenesTable):
            if table.contig not in contig_ids_with_rgp:
                continue
            genedata_to_contig_ids[table.genedata].append(table.contig)
            # contig_genedata_id_to_gene_name[(table.contig, table.genedata)] = table.name

        # Create a contig info dictionary to store contig information
        contig_to_info: dict[int, ContigBorderPosition] = {}

        for table in reader.fetch(GeneDataTable):
            # Problem with RNA genes that are not in GenesTable
            # We could filter based on gene_type columns but need to retrieve this column
            # and convert it to string when parsing GeneDataTable which takes time
            # quick and dirty solution for now:
            try:
                contig_ids = genedata_to_contig_ids[table.idx]
            except KeyError:
                continue
            for contig_id in contig_ids:
                if contig_id not in contig_to_info:
                    contig_to_info[contig_id] = ContigBorderPosition(
                        last_gene_position=table.position,
                        last_gene_idx=table.idx,
                        first_gene_position=table.position,
                        first_gene_idx=table.idx,
                        gene_count=0,
                    )
                contig_to_info[contig_id].gene_count += 1

                if table.position > contig_to_info[contig_id].last_gene_position:
                    contig_to_info[contig_id].last_gene_position = table.position
                    contig_to_info[contig_id].last_gene_idx = table.idx

                if table.position < contig_to_info[contig_id].first_gene_position:
                    contig_to_info[contig_id].first_gene_position = table.position
                    contig_to_info[contig_id].first_gene_idx = table.idx

        assert all(
            info.first_gene_position == 0 for info in contig_to_info.values()
        ), "Some contigs have no gene at position 0"

        contig_id_to_name = {
            contig_info.idx: contig_info.name
            for contig_info in self._rgp_contig_to_info.values()
        }

        # compute mapping (contig_id, genedata_id) to gene_name info only for border genes
        # need to go through GenesTable again to get gene names
        contig_genedata_id_to_gene_name = {
            (contig_id, info.first_gene_idx): None
            for contig_id, info in contig_to_info.items()
        }
        contig_genedata_id_to_gene_name.update(
            {
                (contig_id, info.last_gene_idx): None
                for contig_id, info in contig_to_info.items()
            }
        )
        for table in reader.fetch(GenesTable):
            if (table.contig, table.genedata) in contig_genedata_id_to_gene_name:
                contig_genedata_id_to_gene_name[(table.contig, table.genedata)] = (
                    table.name
                )

        contig_name_to_border_genes = {
            contig_id_to_name[contig_id]: ContigBorderGenes(
                first_gene=contig_genedata_id_to_gene_name[
                    contig_id, info.first_gene_idx
                ],
                last_gene=contig_genedata_id_to_gene_name[
                    contig_id, info.last_gene_idx
                ],
                gene_count=info.gene_count,
            )
            for contig_id, info in contig_to_info.items()
        }

        return contig_name_to_border_genes

    def _get_rgp_genes(self, reader: H5Reader) -> dict[str, set[str]]:
        rgp_genes: dict[str, set[str]] = defaultdict(set)

        for table in reader.fetch(RGPTable):
            rgp_genes[table.rgp].add(table.gene)

        return rgp_genes

    def _get_rgp_info(
        self,
        reader: H5Reader,
        rgp_with_genes: dict[str, set[str]],
    ) -> list[RGPInfo]:

        contig_to_border_genes = self._get_contig_border_genes(reader)

        self._rgp_to_nb_genes: dict[str, int] = {
            rgp_name: len(genes) for rgp_name, genes in rgp_with_genes.items()
        }

        gene_to_family: dict[str, str] = {}
        for table in reader.fetch(GeneFamTable):
            gene_to_family[table.gene] = table.family

        unique_families = set(gene_to_family.values())
        family_ids = {fam: idx for idx, fam in enumerate(unique_families)}

        rgp_infos: list[RGPInfo] = []

        for rgp_name, genes in rgp_with_genes.items():

            rgp_families_ids: BitMap = BitMap()
            rgp_families: set[str] = set()
            contig_name = rgp_name.split("_RGP_")[0]
            for gene in genes:

                fam = gene_to_family[gene]
                fam_id = family_ids[fam]
                rgp_families.add(fam)
                rgp_families_ids.add(fam_id)

            contig_border_info = contig_to_border_genes[contig_name]
            is_contig_circular = self._rgp_contig_to_info[contig_name].is_circular

            is_contig_border = False
            if not is_contig_circular and (
                contig_border_info.first_gene in genes
                or contig_border_info.last_gene in genes
            ):
                is_contig_border = True

            is_whole_contig = False
            if len(genes) == contig_border_info.gene_count:
                is_whole_contig = True

            info = RGPInfo(
                name=rgp_name,
                families=rgp_families,
                families_ids=rgp_families_ids,
                is_contig_border=is_contig_border,
                is_whole_contig=is_whole_contig,
                contig=contig_name,
            )
            rgp_infos.append(info)

        return rgp_infos

    def _construct_single(self, idx: int, rgp: RGPInfo):
        return RegionProxy(
            ID=idx,
            name=rgp.name,
            families=rgp.families_ids,
            modules=BitMap(
                module_id
                for fam in rgp.families
                for module_id in self._fam_to_modules.get(fam, [])
            ),
            contig=rgp.contig,
            organism=self._rgp_contig_to_info[rgp.contig].organism,
            length=self._rgp_to_nb_genes[rgp.name],
            is_contig_border=rgp.is_contig_border,
            is_whole_contig=rgp.is_whole_contig,
        )

    def _construct_single_and_add(self, idx: int, rgp: RGPInfo):
        self.graph.add_node(idx)
        self.rgps.add(self._construct_single(idx, rgp))

    def _construct_multiple(self, idx: int, rgps: list[RGPInfo]):
        return RegionProxy(
            ID=idx,
            name=f"identical_rgps_{self.identical_regions}",
            families=rgps[0].families_ids,
            children=set(
                self._construct_single(i, rgp)
                for i, rgp in enumerate(rgps, start=idx + 1)
            ),
            modules=BitMap(
                module_id
                for fam in rgps[0].families
                for module_id in self._fam_to_modules.get(fam, [])
            ),
            contig=self._contig_to_organism,
            # identical regions object is considered on a contig border if all rgp are contig border
            is_contig_border=all(rgp.is_contig_border for rgp in rgps),
            # identical regions object is considered as whole contig if all rgp are whole contig
            is_whole_contig=all(rgp.is_whole_contig for rgp in rgps),
        )

    def _construct_multiple_and_add(self, idx: int, rgps: list[RGPInfo]):
        self.rgps.add(self._construct_multiple(idx, rgps))
        self.graph.add_node(idx)
        self.identical_regions += 1

    def _construct_and_add(self, idx: int, rgps: list[RGPInfo]):
        if len(rgps) == 1:
            self._construct_single_and_add(idx, rgps[0])
        else:
            self._construct_multiple_and_add(idx, rgps)

    def _grr(self, b1: BitMap, b2: BitMap, mode: Callable) -> float:
        return len(b1 & b2) / mode(len(b1), len(b2))

    def _rgp_metric(self, r1: RegionProxy, r2: RegionProxy, grr_cutoff: float, metric: RGPMetricType) -> RGPMetric:
        if r1.is_contig_border or r2.is_contig_border:
            agrr = self._grr(r1.families, r2.families, min)
            max_grr = self._grr(r1.families, r2.families, max)
            min_grr = agrr
        else:
            agrr = self._grr(r1.families, r2.families, max)
            min_grr = self._grr(r1.families, r2.families, min)
            max_grr = agrr

        m = RGPMetric(max_grr, min_grr, agrr, len(r1.families & r2.families))
        return m if getattr(m, metric) >= grr_cutoff else None

    def _construct_regions(self):
        logging.info("Loading RGPs from pangenome H5 file")

        with H5Reader(self.h5) as reader:

            # self._contig_to_organism = self._get_contig_to_organism(reader)
            # self._contig_id_to_is_circular = self._get_contig_to_is_circular(reader)

            rgp_to_genes = self._get_rgp_genes(reader)

            contigs_with_rgp = {rgp_name.split("_RGP_")[0] for rgp_name in rgp_to_genes}

            self._rgp_contig_to_info = self._get_contig_to_info(
                reader, contigs_with_rgp
            )

            rgp_infos = self._get_rgp_info(reader, rgp_to_genes)

            self._rgp_to_spot = self._get_rgp_spot(reader)
            self._fam_to_modules = self._get_fam_to_modules(reader)

            fams_to_rgps: defaultdict[tuple[int], list[RGPInfo]] = defaultdict(list)

            for info in rgp_infos:
                fams_key = tuple(sorted(fam_id for fam_id in info.families_ids))
                fams_to_rgps[fams_key].append(info)
            print(
                f"{len(fams_to_rgps)} unique family combinations found for {len(rgp_infos)} RGPs"
            )
            idx = 0
            for rgps in fams_to_rgps.values():
                self._construct_and_add(idx, rgps)
                idx += len(rgps) + 1 if len(rgps) > 1 else 1

        logging.info(
            f"{len(rgp_infos)} RGPs loaded from pangenome ({len(self.rgps)} unique RGPs after dereplication)"
        )

    def _compute_all_metrics(self, grr_cutoff: float, metric: RGPMetricType):
        logging.info("Computing RGP metrics")

        nb_pairs = 0
        for r1, r2 in combinations(self.rgps, 2):
            if len(r1.families & r2.families) == 0:
                continue

            nb_pairs += 1
            if m := self._rgp_metric(r1, r2, grr_cutoff, metric):
                self.metrics.append(m)
                self.graph.add_edge(r1.ID, r2.ID, **m.__dict__)
        logging.info(f"RGP metrics computed for {nb_pairs:,} pairs of RGPs ({len(self.metrics)} selected after GRR cutoff)")

    def _louvain_clustering(self, metric: RGPMetricType):
        logging.info(f"Clustering RGPs using Louvain communities on '{metric}' metric")

        partitions = nx.algorithms.community.louvain_communities(
            self.graph, weight=metric
        )
        for i, nodes in enumerate(partitions):
            nx.set_node_attributes(
                self.graph,
                {node: f"cluster_{i}" for node in nodes},
                name=f"{metric}_cluster",
            )

        logging.info(f"Graph has {len(partitions)} clusters using '{metric}'")

    def _add_edges_to_identical_rgps(self):
        logging.info("Unmerging identical RGPs in the graph")

        unmerged = 0
        edge_data = {
            "max_grr": 1.0,
            "min_grr": 1.0,
            "grr": 1.0,
            "identical_famillies": True,
        }

        for rgp in self.rgps:
            if not rgp.is_identical_region:
                continue

            unmerged += 1
            self.graph.add_nodes_from(
                (child.ID for child in rgp.children),
                identical_rcp_group=rgp.name,
            )

            edges = [
                (r1.ID, r2.ID, edge_data)
                for r1, r2 in combinations(rgp.children, 2)
            ]

            for connected in self.graph.neighbors(rgp.ID):
                data = self.graph[rgp.ID][connected]
                edges += [
                    (child.ID, connected, data)
                    for child in rgp.children
                ]

            self.graph.add_edges_from(edges)
            self.graph.remove_node(rgp.ID)

        logging.info(f"Unmerged {unmerged} identical RGPs")

    def _spot_id(self, rgp: RegionProxy) -> str:
        if rgp.name in self._rgp_to_spot:
            return f"spot_{self._rgp_to_spot[rgp.name]}"
        else:
            return "No spot"

    def _add_info_to_identical_rgps(self):
        logging.info("Adding info to identical RGPs in graph")

        identical = 0
        for rgp in self.rgps:
            if not rgp.is_identical_region:
                continue

            identical += 1
            spots = {self._spot_id(child) for child in rgp.children}

            self.graph.nodes[rgp.ID].update({
                "identical_rgp_group": True,
                "name": rgp.name,
                "families_count": len(rgp.families),
                "identical_rgp_count": len(rgp.children),
                "identical_rgp_names": ";".join(child.name for child in rgp.children),
                "identical_rgp_contig_border_count": len(
                    [True for child in rgp.children if child.is_contig_border]
                ),
                "identical_rgp_whole_contig_count": len(
                    [True for child in rgp.children if child.is_whole_contig]
                ),
                "identical_rgp_spots": ";".join(spots),
                "spot_id": (
                    spots.pop()
                    if len(spots) == 1
                    else "Multiple spots"
                ),
                "modules": ",".join(str(module) for module in rgp.modules),
            })

        logging.info(f"Added info to {identical} identical RGPs")

    def _make_info_from_rgp(self, rgp: RegionProxy) -> dict:
        return {
            "contig": rgp.contig,
            "genome": rgp.organism,
            "name": rgp.name,
            "genes_count": rgp.length,
            "is_contig_border": rgp.is_contig_border,
            "is_whole_contig": rgp.is_whole_contig,
            "spot_id": self._spot_id(rgp),
            "modules": ";".join(str(module) for module in rgp.modules),
            "families_count": rgp.nb_families,
        }

    def _add_info_to_rgps(self):
        logging.info("Adding info to RGPs in graph")
        annotated = 0
        for rgp in self.rgps:
            if rgp.ID in self.graph:
                self.graph.nodes[rgp.ID].update(self._make_info_from_rgp(rgp))
                annotated += 1

            if rgp.children:
                for child in rgp.children:
                    if child.ID in self.graph:
                        self.graph.nodes[child.ID].update(self._make_info_from_rgp(child))
                        annotated += 1

        logging.info(f"Added info to {annotated} RGPs")

    def _write_graphs(self, output: Path, basename: str, graph_formats: list[str]):
        if "gexf" in graph_formats:
            graph_filename = output / f"{basename}.gexf"
            logging.info(f"Writing RGP graph in GEXF format to {graph_filename}")
            nx.write_gexf(self.graph, graph_filename)
        if "graphml" in graph_formats:
            graph_filename = output / f"{basename}.graphml"
            logging.info(f"Writing RGP graph in GraphML format to {graph_filename}")
            nx.write_graphml(self.graph, graph_filename)

    def _write_outputs(self, output: Path, basename: str, graph_formats: list[str]):
        self._write_graphs(output, basename, graph_formats)

    def run(self, options: RGPClusteringOptions):
        self._construct_regions()
        self._compute_all_metrics(options.grr_cutoff, options.metric)
        self._louvain_clustering(options.metric)

        if options.unmerge_identical_rgps:
            self._add_edges_to_identical_rgps()
        else:
            self._add_info_to_identical_rgps()

        self._add_info_to_rgps()

        self._write_outputs(options.output, options.basename, options.graph_formats)

    @property
    def rgp_count(self) -> int:
        return len(self.rgps)

class IdenticalRegions:
    """
    Represents a group of Identical Regions within a pangenome.

    :param name: The name of the identical region group.
    :param identical_rgps: A set of Region objects representing the identical regions.
    :param families: A set of GeneFamily objects associated with the identical regions.
    :param is_contig_border: A boolean indicating if the identical regions span across contig borders.
    """

    def __init__(
        self,
        name: str,
        identical_rgps: Set[Region],
        families: Set[GeneFamily],
        is_contig_border: bool,
    ):
        if not isinstance(identical_rgps, set):
            raise TypeError("Expected 'identical_rgps' to be a set")
        else:
            if len(identical_rgps) == 0:
                raise ValueError("Set of identical_rgps must not be empty")
            if not all(isinstance(region, Region) for region in identical_rgps):
                raise TypeError("All element in identical_rgps must be `Region`")
        if not isinstance(families, set):
            raise TypeError("Expected 'families' to be a set")
        else:
            if len(families) == 0:
                raise ValueError("Set of families must not be empty")
            if not all(isinstance(family, GeneFamily) for family in families):
                raise TypeError("All element in families must be `GeneFamilies`")
        self.name = name
        self.families = families
        self.rgps = identical_rgps
        self.is_contig_border = is_contig_border
        self.ID = Region.id_counter

        Region.id_counter += 1

    def __eq__(self, other: "IdenticalRegions") -> bool:
        """
        Check if two IdenticalRegions objects are equal based on their families,
        identical regions, and contig border status.

        :param other: The IdenticalRegions object to compare.
        :return: True if the objects are equal, False otherwise.
        """
        if not isinstance(other, IdenticalRegions):
            # don't attempt to compare against unrelated types
            raise TypeError(
                "'IdenticalRegions' type object was expected, "
                f"but '{type(other)}' type object was provided."
            )

        return (
            self.families == other.families
            and self.rgps == other.rgps
            and self.is_contig_border == other.is_contig_border
        )

    def __repr__(self):
        return (
            f"IdenticalRegions(name='{self.name}', num_rgps={len(self.rgps)}, num_families={len(self.families)},"
            f" is_contig_border={self.is_contig_border})"
        )

    def __str__(self):
        return self.name

    def __hash__(self):
        return id(self)

    def __lt__(self, obj):
        return self.ID < obj.ID

    def __gt__(self, obj):
        return self.ID > obj.ID

    def __le__(self, obj):
        return self.ID <= obj.ID

    def __ge__(self, obj):
        return self.ID >= obj.ID

    @property
    def genes(self):
        """
        Return iterable of genes from all RGPs that are identical in families
        """
        for rgp in self.rgps:
            yield from rgp.genes

    @property
    def spots(self) -> Set[Spot]:
        """
        Return spots from all RGPs that are identical in families
        """
        spots = {rgp.spot for rgp in self.rgps if rgp.spot is not None}
        return spots

    @property
    def modules(self) -> Set[Module]:
        """
        Return iterable of genes from all RGPs that are identical in families
        """
        modules = set()
        for rgp in self.rgps:
            modules |= rgp.modules

        return modules


def compute_grr(
    rgp_a_families: Set[GeneFamily], rgp_b_families: Set[GeneFamily], mode: Callable
) -> float:
    """
    Compute gene repertoire relatedness (GRR) between two rgp.
    Mode can be the function min to compute min GRR or max to compute max_grr

    :param rgp_a_families: Rgp A
    :param rgp_b_families: rgp B
    :param mode: min or max function

    :return: GRR value between 0 and 1
    """

    grr = len(rgp_a_families & rgp_b_families) / mode(
        len(rgp_a_families), len(rgp_b_families)
    )

    return grr


def compute_jaccard_index(rgp_a_families: set, rgp_b_families: set) -> float:
    """
    Compute jaccard index between two rgp based on their families.

    :param rgp_a_families: Rgp A
    :param rgp_b_families: rgp B

    :return : Jaccard index
    """

    jaccard_index = len(rgp_a_families & rgp_b_families) / len(
        rgp_a_families | rgp_b_families
    )

    return jaccard_index


def add_info_to_rgp_nodes(graph, regions: List[Region], region_to_spot: dict):
    """
    Format RGP information into a dictionary for adding to the graph.

    This function takes a list of RGPs and a dictionary mapping each RGP to its corresponding spot ID,
    and formats the RGP information into a dictionary for further processing or addition to a graph.

    :param graph: RGPs graph
    :param regions: A list of RGPs.
    :param region_to_spot: A dictionary mapping each RGP to its corresponding spot ID.
    :return: A dictionary with RGP id as the key and a dictionary containing information on the corresponding RGP as value.
    """

    region_attributes = {}
    for region in regions:
        region_info = {
            "contig": region.contig.name,
            "genome": region.organism.name,
            "name": region.name,
            "genes_count": len(region),
            "is_contig_border": region.is_contig_border,
            "is_whole_contig": region.is_whole_contig,
            "spot_id": get_spot_id(region, region_to_spot),
            "modules": ";".join({str(module) for module in region.modules}),
            "families_count": region.number_of_families,
        }

        region_attributes[region.ID] = region_info

        node_attributes = graph.nodes[region.ID]
        node_attributes.update(region_info)

    return region_attributes


def join_dicts(dicts: List[Dict[str, Any]], delimiter: str = ";") -> Dict[str, Any]:
    """
    Join dictionaries by concatenating the values with a custom delimiter for common keys.

    Given a list of dictionaries, this function creates a new dictionary where the values for common keys
    are concatenated with the specified delimiter.

    :param dicts: A list of dictionaries to be joined.
    :param delimiter: The delimiter to use for joining values. Default is ';'.
    :return: A dictionary with joined values for common keys.
    """
    final_dict = defaultdict(list)
    for dict_obj in dicts:
        for k, v in dict_obj.items():
            final_dict[k].append(str(v))
    return {k: delimiter.join(v) for k, v in final_dict.items()}


def format_rgp_metadata(rgp: Region) -> Dict[str, str]:
    """
    Format RGP metadata by combining source and field values.

    Given an RGP object with metadata, this function creates a new dictionary where the keys
    are formatted as 'source_field' and the values are concatenated with '|' as the delimiter.

    :param rgp: The RGP object with metadata.
    :return: A dictionary with formatted metadata.
    """
    source_field_2_value = defaultdict(list)
    for rgp_metadata in rgp.metadata:
        source = rgp_metadata.source
        for field in rgp_metadata.fields:
            source_field_2_value[f"{source}_{field}"].append(
                str(rgp_metadata.get(field))
            )

    return {
        col_name: "|".join(values) for col_name, values in source_field_2_value.items()
    }


def add_rgp_metadata_to_graph(
    graph: nx.Graph, rgps: List[Union[Region, IdenticalRegions]]
) -> None:
    """
    Add metadata from Region or IdenticalRegions objects to the graph.

    :param graph: The graph to which the metadata will be added.
    :param rgps: A set of Region or IdenticalRegions objects containing the metadata to be added.

    """
    for rgp in rgps:
        element_to_metadata_sources = {
            "family": set(),
            "gene": set(),
            "module": set(),
            "spot": set(),
        }

        for family in rgp.families:
            element_to_metadata_sources["family"] |= {
                metadata.source for metadata in family.metadata
            }
            if family.module:
                element_to_metadata_sources["module"] |= {
                    metadata.source for metadata in family.module.metadata
                }

        for gene in rgp.genes:
            element_to_metadata_sources["gene"] |= {
                metadata.source for metadata in gene.metadata
            }

        if isinstance(rgp, Region):
            rgp_metadata = rgp.formatted_metadata_dict_to_string()
            if rgp.spot is not None:
                element_to_metadata_sources["spot"] = {
                    metadata.source for metadata in rgp.spot.metadata
                }

        elif isinstance(rgp, IdenticalRegions):
            rgp_metadata_dicts = [
                ident_rgp.formatted_metadata_dict_to_string() for ident_rgp in rgp.rgps
            ]
            rgp_metadata = join_dicts(rgp_metadata_dicts)

            element_to_metadata_sources["spot"] |= {
                metadata.source for spot in rgp.spots for metadata in spot.metadata
            }

        else:
            raise TypeError(
                f"Expect Region or IdenticalRegions object, not {type(rgp)}"
            )

        for element, metadata_sources in element_to_metadata_sources.items():
            for source in metadata_sources:
                graph.nodes[rgp.ID][f"has_{element}_with_{source}"] = True

        for metadata_name, value in rgp_metadata.items():
            graph.nodes[rgp.ID][metadata_name] = value


def add_info_to_identical_rgps(
    rgp_graph: nx.Graph,
    identical_rgps_objects: List[IdenticalRegions],
    rgp_to_spot: Dict[Region, int],
):
    """
    Add identical rgps info in the graph as node attributes.

    :params rgp_graph: Graph with rgp id as node and grr value as edges
    :params rgp_to_identical_rgps: dict with uniq RGP as the key and set of identical rgps as value
    """

    for identical_rgp_obj in identical_rgps_objects:
        spots_of_identical_rgp_obj = {
            get_spot_id(i_rgp, rgp_to_spot) for i_rgp in identical_rgp_obj.rgps
        }

        rgp_graph.add_node(
            identical_rgp_obj.ID,
            identical_rgp_group=True,
            name=identical_rgp_obj.name,
            families_count=len(identical_rgp_obj.families),
            identical_rgp_count=len(identical_rgp_obj.rgps),
            identical_rgp_names=";".join(
                i_rgp.name for i_rgp in identical_rgp_obj.rgps
            ),
            identical_rgp_genomes=";".join(
                {i_rgp.organism.name for i_rgp in identical_rgp_obj.rgps}
            ),
            identical_rgp_contig_border_count=len(
                [True for i_rgp in identical_rgp_obj.rgps if i_rgp.is_contig_border]
            ),
            identical_rgp_whole_contig_count=len(
                [True for i_rgp in identical_rgp_obj.rgps if i_rgp.is_whole_contig]
            ),
            identical_rgp_spots=";".join(spots_of_identical_rgp_obj),
            spot_id=(
                spots_of_identical_rgp_obj.pop()
                if len(spots_of_identical_rgp_obj) == 1
                else "Multiple spots"
            ),
            modules=";".join({str(module) for module in identical_rgp_obj.modules}),
        )


def add_edges_to_identical_rgps(
    rgp_graph: nx.Graph, identical_rgps_objects: List[IdenticalRegions]
):
    """
    Replace identical rgp objects by all identical RGPs it contains.

    :param rgp_graph: The RGP graph to add edges to.
    :param identical_rgps_objects: A dictionary mapping RGPs to sets of identical RGPs.
    """

    identical_edge_data = {
        "grr": 1.0,
        "max_grr": 1.0,
        "min_grr": 1.0,
        "identical_famillies": True,
    }

    added_identical_rgps = []

    for identical_rgp_obj in identical_rgps_objects:

        rgp_graph.add_nodes_from(
            [ident_rgp.ID for ident_rgp in identical_rgp_obj.rgps],
            identical_rgp_group=identical_rgp_obj.name,
        )

        # add edge between identical rgp with metrics at one (perfect score)
        edges_to_add = [
            (rgp_a.ID, rgp_b.ID, identical_edge_data)
            for rgp_a, rgp_b in combinations(identical_rgp_obj.rgps, 2)
        ]

        # replicate all edges that connect identical rgp object to other rgps
        for connected_rgp in rgp_graph.neighbors(identical_rgp_obj.ID):
            edge_data = rgp_graph[identical_rgp_obj.ID][connected_rgp]
            edges_to_add += [
                (identical_rgp.ID, connected_rgp, edge_data)
                for identical_rgp in identical_rgp_obj.rgps
            ]

        rgp_graph.add_edges_from(edges_to_add)

        # remove node of the identical rgp object
        rgp_graph.remove_node(identical_rgp_obj.ID)

        added_identical_rgps += list(identical_rgp_obj.rgps)

    return added_identical_rgps


def dereplicate_rgp(
    rgps: Set[Union[Region, IdenticalRegions]], disable_bar: bool = False
) -> List[Union[Region, IdenticalRegions]]:
    """
    Dereplicate RGPs that have the same families.

    Given a list of Region or IdenticalRegions objects representing RGPs, this function groups together
    RGPs with the same families into IdenticalRegions objects and returns a list of dereplicated RGPs.

    :param rgps: A set of Region or IdenticalRegions objects representing the RGPs to be dereplicated.
    :param disable_bar: If True, disable the progress bar.

    :return: A list of dereplicated RGPs (Region or IdenticalRegions objects). For RGPs with the same families,
             they will be grouped together in IdenticalRegions objects.
    """
    logging.info(f"Dereplicating {len(rgps)} RGPs")
    families_to_rgps = defaultdict(list)

    for rgp in tqdm(rgps, total=len(rgps), unit="RGP", disable=disable_bar):
        families_to_rgps[tuple(sorted(f.ID for f in rgp.families))].append(rgp)

    dereplicated_rgps = []
    identical_region_count = 0
    for rgps in families_to_rgps.values():
        if len(rgps) == 1:
            dereplicated_rgps.append(rgps[0])
        else:
            families = set(rgps[0].families)

            # identical regions object is considered on a contig border if all rgp are contig border
            is_contig_border = all(rgp.is_contig_border for rgp in rgps)

            # create a new object that will represent the identical rgps
            identical_rgp = IdenticalRegions(
                name=f"identical_rgps_{identical_region_count}",
                identical_rgps=set(rgps),
                families=families,
                is_contig_border=is_contig_border,
            )
            identical_region_count += 1
            dereplicated_rgps.append(identical_rgp)

    logging.info(f"{len(dereplicated_rgps)} unique RGPs")
    return dereplicated_rgps


def compute_rgp_metric(
    rgp_a: Region, rgp_b: Region, grr_cutoff: float, grr_metric: str
) -> Union[Tuple[int, int, dict], None]:
    """
    Compute GRR metric between two RGPs.

    :param rgp_a: A rgp
    :param rgp_b: another rgp
    :param grr_cutoff: Cutoff filter
    :param grr_metric: grr mode between min_grr, max_grr and incomplete_aware_grr

    :returns: Tuple containing the IDs of the two RGPs and the computed metrics as a dictionary
    """

    edge_metrics = {}

    # RGP at a contig border are seen as incomplete and min GRR is used instead of max GRR
    if rgp_a.is_contig_border or rgp_b.is_contig_border:
        edge_metrics["incomplete_aware_grr"] = compute_grr(
            set(rgp_a.families), set(rgp_b.families), min
        )
    else:
        edge_metrics["incomplete_aware_grr"] = compute_grr(
            set(rgp_a.families), set(rgp_b.families), max
        )

    # Compute max and min GRR metrics
    edge_metrics["max_grr"] = compute_grr(set(rgp_a.families), set(rgp_b.families), max)
    edge_metrics["min_grr"] = compute_grr(set(rgp_a.families), set(rgp_b.families), min)

    # The number of shared families can be useful when visualizing the graph
    edge_metrics["shared_family"] = len(
        set(rgp_a.families).intersection(set(rgp_b.families))
    )

    # Only return the metrics if the GRR value is above the cutoff
    if edge_metrics[grr_metric] >= grr_cutoff:
        return rgp_a.ID, rgp_b.ID, edge_metrics

def cluster_rgp_on_grr(graph: nx.Graph, clustering_attribute: str = "grr"):
    """
    Cluster rgp based on grr using louvain communities clustering.

    :param graph: NetworkX graph object representing the RGPs and their relationship
    :param clustering_attribute: Attribute of the graph to use for clustering (default is "grr")
    """

    partitions = nx.algorithms.community.louvain_communities(
        graph, weight=clustering_attribute
    )

    # Add partition index in node attributes
    for i, cluster_nodes in enumerate(partitions):
        nx.set_node_attributes(
            graph,
            {node: f"cluster_{i}" for node in cluster_nodes},
            name=f"{clustering_attribute}_cluster",
        )

    logging.info(f"Graph has {len(partitions)} clusters using {clustering_attribute}")


def get_spot_id(rgp: Region, rgp_to_spot: Dict[Region, int]) -> str:
    """
    Return Spot ID associated to an RGP.
    It adds the prefix "spot" to the spot ID. When no spot is associated with the RGP,
    then the string "No spot" is return

    :param rgp: RGP id
    :param rgp_to_spot: A dictionary mapping an RGP to its spot.

    :return: Spot ID of the given RGP with the prefix spot or "No spot".
    """
    if rgp in rgp_to_spot:
        return f"spot_{rgp_to_spot[rgp]}"
    else:
        return "No spot"


def write_rgp_cluster_table(
    outfile: str,
    grr_graph: nx.Graph,
    rgps_in_graph: List[Union[Region, IdenticalRegions]],
    grr_metric: str,
    rgp_to_spot: Dict[Region, int],
) -> None:
    """
    Writes RGP cluster info to a TSV file using pandas.

    :param outfile: Name of the tsv file
    :param grr_graph: The GRR graph.
    :param rgps_in_graph: A dictionary mapping an RGP to a set of identical RGPs.
    :param grr_metric: The GRR metric used for clustering.
    :param rgp_to_spot: A dictionary mapping an RGP to its spot.
    :return: None
    """

    all_rgps_infos = []
    for rgp_in_graph in rgps_in_graph:
        cluster = grr_graph.nodes[rgp_in_graph.ID][f"{grr_metric}_cluster"]

        identical_rgps = (
            [rgp_in_graph] if isinstance(rgp_in_graph, Region) else rgp_in_graph.rgps
        )

        all_rgps_infos += [
            {"RGPs": r.name, "cluster": cluster, "spot_id": get_spot_id(r, rgp_to_spot)}
            for r in identical_rgps
        ]

    df = pd.DataFrame(all_rgps_infos)
    df.to_csv(outfile, sep="\t", index=False)


def cluster_rgp(
    pangenome,
    grr_cutoff: float,
    output: str,
    basename: str,
    ignore_incomplete_rgp: bool,
    unmerge_identical_rgps: bool,
    grr_metric: str,
    disable_bar: bool,
    graph_formats: Set[str],
    add_metadata: bool = False,
    metadata_sep: str = "|",
    metadata_sources: List[str] = None,
):
    """
    Main function to cluster regions of genomic plasticity based on their GRR

    :param pangenome: pangenome object
    :param grr_cutoff: GRR cutoff value for clustering
    :param output: Directory where the output files will be saved
    :param basename: Basename for the output files
    :param ignore_incomplete_rgp: Whether to ignore incomplete RGPs located at a contig border
    :param unmerge_identical_rgps: Whether to unmerge identical RGPs into separate nodes in the graph
    :param grr_metric: GRR metric to use for clustering
    :param disable_bar: Whether to disable the progress bar
    :param graph_formats: Set of graph file formats to save the output
    :param add_metadata: Add metadata to cluster files
    :param metadata_sep: The separator used to join multiple metadata values
    :param metadata_sources: Sources of the metadata to use and write in the outputs. None means all sources are used.
    """

    metatypes = set()
    need_metadata = False
    if add_metadata:
        for element in ["RGPs", "genes", "spots", "families", "modules"]:
            if pangenome.status["metadata"][element] == "inFile":

                sources_to_use = set(pangenome.status["metasources"][element])

                if metadata_sources is not None:
                    if (
                        len(
                            set(pangenome.status["metasources"][element])
                            & set(metadata_sources)
                        )
                        == 0
                    ):
                        logging.info(
                            f"Metadata for {element} found in pangenome, but none match the specified sources {metadata_sources}. "
                            f"Current source for {element}: {sources_to_use}."
                        )
                        continue
                    else:
                        sources_to_use = set(
                            pangenome.status["metasources"][element]
                        ) & set(metadata_sources)

                need_metadata = True
                metatypes.add(element)
                logging.info(
                    f"Metadata for {element} found in pangenome with sources {sources_to_use}. They will be included in the RGP graph."
                )

    # check statuses and load info
    check_pangenome_info(
        pangenome,
        need_families=True,
        need_annotations=True,
        disable_bar=disable_bar,
        need_rgp=True,
        need_spots=True,
        need_modules=True,
        need_metadata=need_metadata,
        sources=metadata_sources,
        metatypes=metatypes,
    )

    if pangenome.regions == 0:
        raise Exception(
            "The pangenome has no RGPs. The clustering of RGP is then not possible."
        )

    # add all rgp as node
    if ignore_incomplete_rgp:
        valid_rgps = [rgp for rgp in pangenome.regions if not rgp.is_contig_border]

        ignored_rgp_count = pangenome.number_of_rgp - len(valid_rgps)
        total_rgp_count = pangenome.number_of_rgp

        logging.info(
            f"Ignoring {ignored_rgp_count}/{total_rgp_count} ({100 * ignored_rgp_count / total_rgp_count:.2f}%) "
            "RGPs that are located at a contig border and are likely incomplete."
        )

        if len(valid_rgps) == 0:
            raise Exception(
                "The pangenome has no complete RGPs. The clustering of RGP is then not possible."
            )
    else:
        valid_rgps = set(pangenome.regions)

    dereplicated_rgps = dereplicate_rgp(valid_rgps, disable_bar=disable_bar)

    grr_graph = nx.Graph()
    grr_graph.add_nodes_from(rgp.ID for rgp in dereplicated_rgps)

    # Get all pairs of RGP that share at least one family

    family2rgp = defaultdict(set)
    for rgp in dereplicated_rgps:
        for fam in rgp.families:
            family2rgp[fam].add(rgp)

    rgp_pairs = set()
    for rgps in family2rgp.values():
        rgp_pairs |= {tuple(sorted(rgp_pair)) for rgp_pair in combinations(rgps, 2)}

    pairs_count = len(rgp_pairs)

    logging.info(f"Computing GRR metric for {pairs_count:,} pairs of RGP.")

    pairs_of_rgps_metrics = []

    for rgp_a, rgp_b in rgp_pairs:

        pair_metrics = compute_rgp_metric(rgp_a, rgp_b, grr_cutoff, grr_metric)

        if pair_metrics:
            pairs_of_rgps_metrics.append(pair_metrics)

    grr_graph.add_edges_from(pairs_of_rgps_metrics)
    identical_rgps_objects = [
        rgp for rgp in dereplicated_rgps if isinstance(rgp, IdenticalRegions)
    ]
    rgp_objects_in_graph = [rgp for rgp in dereplicated_rgps if isinstance(rgp, Region)]

    if unmerge_identical_rgps:
        rgp_objects_in_graph += add_edges_to_identical_rgps(
            grr_graph, identical_rgps_objects
        )

    # cluster rgp based on grr value
    logging.info(
        f"Louvain_communities clustering of RGP  based on {grr_metric} on {grr_graph}."
    )

    cluster_rgp_on_grr(grr_graph, grr_metric)

    rgp_to_spot = {
        region: int(spot.ID) for spot in pangenome.spots for region in spot.regions
    }

    if not unmerge_identical_rgps:
        logging.info("Add info on identical RGPs merged in the graph")
        add_info_to_identical_rgps(grr_graph, identical_rgps_objects, rgp_to_spot)

    rgps_in_graph = (
        rgp_objects_in_graph if unmerge_identical_rgps else dereplicated_rgps
    )

    # add some attribute to the graph nodes.
    logging.info("Add RGP information to the graph")
    add_info_to_rgp_nodes(grr_graph, rgp_objects_in_graph, rgp_to_spot)

    if need_metadata:
        add_rgp_metadata_to_graph(grr_graph, rgps_in_graph)

    if "gexf" in graph_formats:
        # writing graph in gexf format
        graph_file_name = os.path.join(output, f"{basename}.gexf")
        logging.info(f"Writing graph in gexf format in {graph_file_name}.")
        nx.readwrite.gexf.write_gexf(grr_graph, graph_file_name)

    if "graphml" in graph_formats:
        graph_file_name = os.path.join(output, f"{basename}.graphml")
        logging.info(f"Writing graph in graphml format in {graph_file_name}.")
        nx.readwrite.graphml.write_graphml(grr_graph, graph_file_name)

    outfile = os.path.join(output, f"{basename}.tsv")
    logging.info(f"Writing rgp clusters in tsv format in {outfile}")

    write_rgp_cluster_table(outfile, grr_graph, rgps_in_graph, grr_metric, rgp_to_spot)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provided by user
    """
    pangenome = Pangenome()

    mk_outdir(args.output, args.force)

    pangenome.add_file(args.pangenome)

    A = int(os.environ.get("PPNEW"))

    if A == 0 or A == 2:
        with Timer("OLD", logging):
            cluster_rgp(
                pangenome,
                grr_cutoff=args.grr_cutoff,
                output=args.output,
                basename=args.basename,
                ignore_incomplete_rgp=args.ignore_incomplete_rgp,
                unmerge_identical_rgps=args.no_identical_rgp_merging,
                grr_metric=args.grr_metric,
                disable_bar=args.disable_prog_bar,
                graph_formats=args.graph_formats,
                add_metadata=args.add_metadata,
                metadata_sep=args.metadata_sep,
                metadata_sources=args.metadata_sources,
            )
    if A == 1 or A == 2: 
        with Timer("NEW", logging):
            clustering = RGPClustering(args.pangenome)
            clustering.run(
                RGPClusteringOptions(
                    unmerge_identical_rgps=args.no_identical_rgp_merging,
                    grr_cutoff=args.grr_cutoff,
                    metric=args.grr_metric,
                    output=args.output,
                    basename=args.basename,
                    graph_formats=args.graph_formats,
                    with_metadata=args.add_metadata,
                    metadata_sep=args.metadata_sep,
                    metadata_sources=args.metadata_sources,
                )
            )
    if A == 2:
        # comparison of RGP between old and new implementation
        rgp_name_to_rgp_proxy = {}
        for rgp_proxy in clustering.rgps:
            if rgp_proxy.children:
                for child_rgp_proxy in rgp_proxy.children:
                    rgp_name_to_rgp_proxy[child_rgp_proxy.name] = child_rgp_proxy
            else:
                rgp_name_to_rgp_proxy[rgp_proxy.name] = rgp_proxy

        for region in pangenome.regions:
            rgp_proxy = rgp_name_to_rgp_proxy[region.name]

            # Log debug info when either region has is_contig_border=True or is_whole_contig=True
            if region.is_contig_border or rgp_proxy.is_contig_border or region.is_whole_contig or rgp_proxy.is_whole_contig:
                logging.debug(
                    f"Comparing RGP: {region.name}\n"
                    f"  Region object:\n"
                    f"    - name: {region.name}\n"
                    f"    - is_contig_border: {region.is_contig_border}\n"
                    f"    - is_whole_contig: {region.is_whole_contig}\n"
                    f"    - families count: {region.number_of_families}\n"
                    f"    - genes count: {len(region)}\n"
                    f"    - contig: {region.contig.name}\n"
                    f"    - organism: {region.organism.name}\n"
                    f"  RegionProxy object:\n"
                    f"    - name: {rgp_proxy.name}\n"
                    f"    - is_contig_border: {rgp_proxy.is_contig_border}\n"
                    f"    - is_whole_contig: {rgp_proxy.is_whole_contig}\n"
                    f"    - families count: {len(rgp_proxy.families)}\n"
                    f"    - genes count: {rgp_proxy.length}\n"
                    f"    - contig: {rgp_proxy.contig}\n"
                    f"    - organism: {rgp_proxy.organism}"
                )

            if region.is_contig_border != rgp_proxy.is_contig_border:
                logging.error(f"Mismatch in is_contig_border for RGP: {region.name}")

            if region.is_whole_contig != rgp_proxy.is_whole_contig:
                logging.error(f"Mismatch in is_whole_contig for RGP: {region.name}")

            assert region.is_contig_border == rgp_proxy.is_contig_border, region.name
            assert region.is_whole_contig == rgp_proxy.is_whole_contig, region.name


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : Sub_parser for cluster_rgp command

    :return : Parser arguments for cluster_rgp command
    """
    parser = sub_parser.add_parser(
        "rgp_cluster", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_cluster_rgp(parser)
    return parser


def parser_cluster_rgp(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of rgp command

    :param parser: Parser for cluster_rgp argument
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "-p", "--pangenome", required=True, type=Path, help="The pangenome .h5 file"
    )

    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument(
        "--grr_cutoff",
        required=False,
        type=restricted_float,
        default=0.8,
        help="Min gene repertoire relatedness metric used in the rgp clustering",
    )
    optional.add_argument(
        "--grr_metric",
        required=False,
        type=str,
        default="incomplete_aware_grr",
        help="The grr (Gene Repertoire Relatedness) is used to assess the similarity between two "
        "RGPs based on their gene families."
        "There are three different modes for calculating the grr value: 'min_grr', 'max_grr' "
        "or  'incomplete_aware_grr'."
        "'min_grr': Computes the number of gene families shared between the two RGPs and "
        "divides it by the smaller number of gene families among the two RGPs."
        "'max_grr': Calculates the number of gene families shared between the two RGPs and "
        "divides it by the larger number of gene families among the two RGPs."
        "'incomplete_aware_grr' (default): If at least one RGP is considered incomplete, "
        "which occurs when it is located at the border of a contig,"
        "the 'min_grr' mode is used. Otherwise, the 'max_grr' mode is applied.",
        choices=["incomplete_aware_grr", "min_grr", "max_grr"],
    )

    optional.add_argument(
        "--ignore_incomplete_rgp",
        required=False,
        action="store_true",
        help="Do not cluster RGPs located on a contig border which are likely incomplete.",
    )

    optional.add_argument(
        "--no_identical_rgp_merging",
        required=False,
        action="store_true",
        help="Do not merge in one node identical RGP "
        "(i.e. having the same family content) before clustering.",
    )

    optional.add_argument(
        "--basename",
        required=False,
        default="rgp_cluster",
        help="basename for the output file",
    )

    optional.add_argument(
        "-o",
        "--output",
        required=False,
        type=Path,
        default="rgp_clustering",
        help="Output directory",
    )

    optional.add_argument(
        "--graph_formats",
        required=False,
        type=str,
        choices=["gexf", "graphml"],
        nargs="+",
        default=["gexf", "graphml"],
        help="Format of the output graph.",
    )

    optional.add_argument(
        "--add_metadata",
        required=False,
        action="store_true",
        help="Include metadata information in the output files "
        "if any have been added to pangenome elements (see ppanggolin metadata command).",
    )

    optional.add_argument(
        "--metadata_sources",
        default=None,
        nargs="+",
        help="Which source of metadata should be written. "
        "By default all metadata sources are included.",
    )

    optional.add_argument(
        "--metadata_sep",
        required=False,
        default="|",
        help="The separator used to join multiple metadata values for elements with multiple metadata"
        " values from the same source. This character should not appear in metadata values.",
    )
