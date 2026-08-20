#!/usr/bin/env python3

# default libraries
import logging
import argparse
import json
from itertools import combinations
from collections.abc import Callable
from collections import defaultdict
from typing import Dict, List, Set, Union
from pathlib import Path

# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd

# local libraries
from ppanggolin.region import Region, Spot, Module
from ppanggolin.utils import restricted_float, mk_outdir
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
        self.modules: set[int] = set(modules) if modules is not None else set()

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
        # Hash on the stable, deterministically-assigned ID rather than id(self):
        # object addresses vary between runs and make set/dict iteration order
        # (and therefore graph construction order) non-reproducible.
        return hash(self.ID)

    def __eq__(self, rhs: "RegionProxy") -> bool:
        return (
            self.families == rhs.families
            and self.children == rhs.children
            and self.is_contig_border == rhs.is_contig_border
        )

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
    rgp: tp.Annotated[
        str, TableAttribute(name="RGP", transform=lambda x: x.decode("utf-8"))
    ]
    gene: tp.Annotated[
        str, TableAttribute(name="gene", transform=lambda x: x.decode("utf-8"))
    ]
    _table: str = "/RGP"


@dataclass
class GeneFamTable:
    gene: tp.Annotated[
        str, TableAttribute(name="gene", transform=lambda x: x.decode("utf-8"))
    ]
    family: tp.Annotated[
        str, TableAttribute(name="geneFam", transform=lambda x: x.decode("utf-8"))
    ]
    _table: str = "/geneFamilies"


@dataclass
class AnnotationsGeneTable:
    name: tp.Annotated[
        str, TableAttribute(name="ID", transform=lambda x: x.decode("utf-8"))
    ]
    contig: tp.Annotated[int, TableAttribute(name="contig")]
    _table: str = "/annotations/genes"


@dataclass
class RGPSpotTable:
    rgp: tp.Annotated[
        str, TableAttribute(name="RGP", transform=lambda x: x.decode("utf-8"))
    ]
    spot: tp.Annotated[int, TableAttribute(name="spot")]
    _table: str = "/spots"


@dataclass
class ModuleTable:
    fam: tp.Annotated[
        str, TableAttribute(name="geneFam", transform=lambda x: x.decode("utf-8"))
    ]
    module: tp.Annotated[int, TableAttribute(name="module")]
    _table: str = "/modules"


@dataclass
class ContigTable:
    genome: tp.Annotated[
        str, TableAttribute(name="genome", transform=lambda x: x.decode("utf-8"))
    ]
    contig: tp.Annotated[
        str, TableAttribute(name="name", transform=lambda x: x.decode("utf-8"))
    ]
    is_circular: tp.Annotated[bool, TableAttribute(name="is_circular")]
    idx: tp.Annotated[int, TableAttribute(name="ID")]
    _table: str = "/annotations/contigs"


@dataclass
class GenesTable:
    name: tp.Annotated[
        str, TableAttribute(name="ID", transform=lambda x: x.decode("utf-8"))
    ]
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
    ignore_incomplete_rgp: bool = False
    disable_prog_bar: bool = False


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
        self._rgp_to_genes: dict[str, set[str]] = None
        self._family_id_to_name: dict[int, str] = None
        self._has_metadata: bool = False
        self._rgp_metadata_values: Dict[str, Dict[str, list]] = defaultdict(
            lambda: defaultdict(list)
        )
        self._family_metadata_sources: Dict[str, set] = defaultdict(set)
        self._gene_metadata_sources: Dict[str, set] = defaultdict(set)
        self._spot_metadata_sources: Dict[int, set] = defaultdict(set)
        self._module_metadata_sources: Dict[int, set] = defaultdict(set)

    def _check_required_status(
        self,
        status: dict,
        required_statuses: dict[str, str],
    ):
        """Check that all required pangenome information is available in the H5 file."""

        for status_name, description in required_statuses.items():
            if status.get(status_name) != "inFile":
                raise ValueError(
                    f"Cannot compute RGP metrics: {description} is not available "
                    "in the pangenome H5 file."
                )

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

    def _fetch_required_table(self, reader: H5Reader, table_type):
        try:
            return reader.fetch(table_type)
        except KeyError as e:
            raise ValueError(
                f"Required table '{table_type.__name__}' is missing from the "
                "pangenome H5 file."
            ) from e

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

    def _get_metadata_sources(
        self, reader: H5Reader, metadata_sources: list[str] = None
    ) -> dict[str, list[str]]:
        """
        Get, per metatype, the list of metadata sources present in the pangenome file,
        restricted to the requested sources (if any).
        """
        status_attrs = reader.handle.root.status._v_attrs
        if not (hasattr(status_attrs, "metadata") and status_attrs.metadata):
            return {}

        metastatus = reader.handle.root.status.metastatus._v_attrs
        metasources = reader.handle.root.status.metasources._v_attrs

        metatype_to_sources: dict[str, list[str]] = {}
        for metatype in ("RGPs", "genes", "spots", "families", "modules"):
            if not getattr(metastatus, metatype, False):
                continue

            sources = list(getattr(metasources, metatype, []))
            if metadata_sources is not None:
                sources = [source for source in sources if source in metadata_sources]

            if sources:
                metatype_to_sources[metatype] = sources
                logging.getLogger("PPanGGOLiN").info(
                    f"Metadata for {metatype} found in pangenome with sources {sources}. "
                    "They will be included in the RGP graph."
                )
            elif metadata_sources is not None:
                logging.getLogger("PPanGGOLiN").info(
                    f"Metadata for {metatype} found in pangenome, but none match the "
                    f"specified sources {metadata_sources}."
                )

        return metatype_to_sources

    def _read_metadata_source(
        self,
        reader: H5Reader,
        metatype: str,
        source: str,
        value_target: Dict[str, Dict[str, list]] = None,
        source_target: Dict[Union[str, int], set] = None,
    ) -> None:
        """
        Read a single metadata source table directly from the H5 file (lightweight, no
        pangenome object creation), and accumulate either the raw values (for RGPs metadata,
        which are written as node attributes) and/or the set of sources providing metadata
        for each identifier (used for family/gene/module/spot presence flags).
        """
        table = reader.handle.get_node(f"/metadata/{metatype}/{source}")

        for row in table.iterrows():
            identifier = row["ID"]
            if isinstance(identifier, bytes):
                identifier = identifier.decode("utf-8")

            if source_target is not None:
                source_target[identifier].add(source)

            if value_target is not None:
                for field in table.colnames:
                    if field == "ID":
                        continue
                    value = row[field]
                    if isinstance(value, bytes):
                        value = value.decode("utf-8")
                    value_target[identifier][f"{source}_{field}"].append(str(value))

    def _load_metadata(
        self, reader: H5Reader, metadata_sources: list[str] = None
    ) -> None:
        """
        Load metadata directly from the H5 file, in the same lightweight way RGP info is
        loaded (reading tables directly instead of building full pangenome objects).
        """
        metatype_to_sources = self._get_metadata_sources(reader, metadata_sources)

        for metatype, sources in metatype_to_sources.items():
            for source in sources:
                if metatype == "RGPs":
                    self._read_metadata_source(
                        reader, metatype, source, value_target=self._rgp_metadata_values
                    )
                elif metatype == "families":
                    self._read_metadata_source(
                        reader,
                        metatype,
                        source,
                        source_target=self._family_metadata_sources,
                    )
                elif metatype == "genes":
                    self._read_metadata_source(
                        reader,
                        metatype,
                        source,
                        source_target=self._gene_metadata_sources,
                    )
                elif metatype == "spots":
                    self._read_metadata_source(
                        reader,
                        metatype,
                        source,
                        source_target=self._spot_metadata_sources,
                    )
                elif metatype == "modules":
                    self._read_metadata_source(
                        reader,
                        metatype,
                        source,
                        source_target=self._module_metadata_sources,
                    )

        self._has_metadata = bool(metatype_to_sources)

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

        # Sort before enumerating: iterating a set of strings depends on
        # PYTHONHASHSEED, which is randomized per process by default and would
        # make family ID assignment (and therefore node ordering) non-reproducible.
        unique_families = sorted(set(gene_to_family.values()))
        family_ids = {fam: idx for idx, fam in enumerate(unique_families)}
        self._family_id_to_name = {idx: fam for fam, idx in family_ids.items()}

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

    def _rgp_metric(
        self, r1: RegionProxy, r2: RegionProxy, grr_cutoff: float, metric: RGPMetricType
    ) -> RGPMetric:
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

    def _construct_regions(
        self,
        with_metadata: bool = False,
        metadata_sources: list[str] = None,
        ignore_incomplete_rgp: bool = False,
    ):
        logging.getLogger("PPanGGOLiN").info("Loading RGPs from pangenome H5 file")

        with H5Reader(self.h5) as reader:

            rgp_to_genes = self._get_rgp_genes(reader)
            self._rgp_to_genes = rgp_to_genes

            contigs_with_rgp = {rgp_name.split("_RGP_")[0] for rgp_name in rgp_to_genes}

            self._rgp_contig_to_info = self._get_contig_to_info(
                reader, contigs_with_rgp
            )

            rgp_infos = self._get_rgp_info(reader, rgp_to_genes)

            if ignore_incomplete_rgp:
                rgp_infos = [
                    rgp_info for rgp_info in rgp_infos if not rgp_info.is_contig_border
                ]
                logging.getLogger("PPanGGOLiN").info(
                    f"{len(rgp_infos)} RGPs loaded from pangenome after filtering out incomplete RGPs"
                )

            self._rgp_to_spot = self._get_rgp_spot(reader)
            self._fam_to_modules = self._get_fam_to_modules(reader)

            if with_metadata:
                self._load_metadata(reader, metadata_sources)

            fams_to_rgps: defaultdict[tuple[int], list[RGPInfo]] = defaultdict(list)

            for info in rgp_infos:
                fams_key = tuple(sorted(fam_id for fam_id in info.families_ids))
                fams_to_rgps[fams_key].append(info)

            idx = 0
            for _fams, rgps in sorted(fams_to_rgps.items(), key=lambda x: x[0]):
                self._construct_and_add(idx, rgps)
                idx += len(rgps) + 1 if len(rgps) > 1 else 1

        logging.getLogger("PPanGGOLiN").info(
            f"{len(rgp_infos)} RGPs loaded from pangenome ({len(self.rgps)} unique RGPs after dereplication)"
        )
        if len(rgp_infos) == 0:
            logging.getLogger("PPanGGOLiN").warning("Pangenome contains no RGPs. Output files will be empty.")

    def _compute_all_metrics(
        self, grr_cutoff: float, metric: RGPMetricType, disable_bar: bool = False
    ):
        """Compute metrics without materializing the set of candidate pairs."""
        logging.getLogger("PPanGGOLiN").info("Computing RGP metrics")

        family_to_rgps = defaultdict(list)
        for rgp in self.rgps:
            for family in rgp.families:
                family_to_rgps[family].append(rgp)

        nb_pairs = 0

        with tqdm(
            total=len(family_to_rgps),
            desc="Computing metrics",
            disable=disable_bar,
            bar_format="{desc}: {percentage:3.0f}%|{bar}",
        ) as pbar:
            for family in sorted(family_to_rgps):
                for r1, r2 in combinations(family_to_rgps[family], 2):
                    shared_families = r1.families & r2.families
                    if family != next(iter(shared_families)):
                        continue

                    nb_pairs += 1
                    if m := self._rgp_metric(r1, r2, grr_cutoff, metric):
                        self.metrics.append(m)
                        self.graph.add_edge(r1.ID, r2.ID, **m.__dict__)
                pbar.update(1)

        logging.getLogger("PPanGGOLiN").info(
            f"RGP metrics computed for {nb_pairs:,} pairs of RGPs "
            f"({len(self.metrics):,} selected after GRR cutoff)"
        )


    def _louvain_clustering(self, metric: RGPMetricType):
        logging.getLogger("PPanGGOLiN").info(f"Clustering RGPs using Louvain communities on '{metric}' metric")

        partitions = nx.algorithms.community.louvain_communities(
            self.graph, weight=metric, seed=42
        )
        ordered_partitions = sorted(
            partitions,
            key=lambda nodes: (min(nodes), len(nodes)),
        )

        for i, nodes in enumerate(ordered_partitions):
            nx.set_node_attributes(
                self.graph,
                {node: f"cluster_{i}" for node in nodes},
                name=f"{metric}_cluster",
            )

        logging.getLogger("PPanGGOLiN").info(f"Graph has {len(ordered_partitions)} RGP clusters using '{metric}'")

    def _add_edges_to_identical_rgps(self):
        logging.getLogger("PPanGGOLiN").info("Unmerging identical RGPs in the graph")

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
                identical_rgp_group=rgp.name,
            )

            edges = [
                (r1.ID, r2.ID, edge_data) for r1, r2 in combinations(rgp.children, 2)
            ]

            for connected in self.graph.neighbors(rgp.ID):
                data = self.graph[rgp.ID][connected]
                edges += [(child.ID, connected, data) for child in rgp.children]

            self.graph.add_edges_from(edges)
            self.graph.remove_node(rgp.ID)

        logging.getLogger("PPanGGOLiN").info(f"Unmerged {unmerged} identical RGPs")

    def _spot_id(self, rgp: RegionProxy) -> str:
        if rgp.name in self._rgp_to_spot:
            return f"spot_{self._rgp_to_spot[rgp.name]}"
        else:
            return "No spot"

    def _add_metadata_info(
        self,
        info: dict,
        rgp_names: List[str],
        families_ids: BitMap,
        module_ids: Set[int],
    ) -> None:
        """
        Add metadata-derived attributes to a node info dict: booleans indicating whether
        family/gene/module/spot metadata is present from a given source, and the RGP's own
        metadata fields.
        """
        if not self._has_metadata:
            return

        element_to_sources = {
            "family": (
                set().union(
                    *(
                        self._family_metadata_sources.get(
                            self._family_id_to_name[fid], set()
                        )
                        for fid in families_ids
                    )
                )
                if families_ids
                else set()
            ),
            "module": (
                set().union(
                    *(
                        self._module_metadata_sources.get(module_id, set())
                        for module_id in module_ids
                    )
                )
                if module_ids
                else set()
            ),
            "gene": set(),
            "spot": set(),
        }

        for rgp_name in rgp_names:
            for gene_name in self._rgp_to_genes.get(rgp_name, ()):
                element_to_sources["gene"] |= self._gene_metadata_sources.get(
                    gene_name, set()
                )

            spot_id = self._rgp_to_spot.get(rgp_name)
            if spot_id is not None:
                element_to_sources["spot"] |= self._spot_metadata_sources.get(
                    spot_id, set()
                )

        for element, sources in element_to_sources.items():
            for source in sources:
                info[f"has_{element}_with_{source}"] = True

        combined_rgp_metadata: Dict[str, list] = defaultdict(list)
        for rgp_name in rgp_names:
            for key, values in self._rgp_metadata_values.get(rgp_name, {}).items():
                combined_rgp_metadata[key].extend(values)

        for key, values in combined_rgp_metadata.items():
            info[key] = json.dumps(values)

    def _make_info_identical_rgp(self, rgp: RegionProxy) -> dict:

        spots = {self._spot_id(child) for child in rgp.children}

        info = {
            "identical_rgp_group": True,
            "name": rgp.name,
            "families_count": len(rgp.families),
            "identical_rgp_count": len(rgp.children),
            "identical_rgp_names": ";".join(child.name for child in rgp.children),
            "identical_rgp_genomes": ";".join(
                {child.organism for child in rgp.children}
            ),
            "identical_rgp_contig_border_count": len(
                [True for child in rgp.children if child.is_contig_border]
            ),
            "identical_rgp_whole_contig_count": len(
                [True for child in rgp.children if child.is_whole_contig]
            ),
            "identical_rgp_spots": ";".join(spots),
            "spot_id": (spots.pop() if len(spots) == 1 else "Multiple spots"),
            "modules": ";".join(f"module_{module}" for module in rgp.modules),
        }

        self._add_metadata_info(
            info,
            [child.name for child in rgp.children],
            rgp.families,
            rgp.modules,
        )

        return info

    def _make_info_from_rgp(self, rgp: RegionProxy) -> dict:

        info = {
            "contig": rgp.contig,
            "genome": rgp.organism,
            "name": rgp.name,
            "genes_count": rgp.length,
            "is_contig_border": rgp.is_contig_border,
            "is_whole_contig": rgp.is_whole_contig,
            "spot_id": self._spot_id(rgp),
            "modules": ";".join(f"module_{module}" for module in rgp.modules),
            "families_count": rgp.nb_families,
        }

        self._add_metadata_info(info, [rgp.name], rgp.families, rgp.modules)

        return info

    def _add_info_to_rgps(self):

        logging.getLogger("PPanGGOLiN").info("Adding info to RGPs in graph")

        annotated = 0
        for rgp in self.rgps:

            if rgp.ID in self.graph:
                if rgp.is_identical_region:
                    self.graph.nodes[rgp.ID].update(self._make_info_identical_rgp(rgp))

                else:
                    self.graph.nodes[rgp.ID].update(self._make_info_from_rgp(rgp))

                annotated += 1

            if rgp.children:  # in case identical rgp are unmerged
                for child in rgp.children:
                    if child.ID in self.graph:
                        self.graph.nodes[child.ID].update(
                            self._make_info_from_rgp(child)
                        )
                        annotated += 1

        logging.getLogger("PPanGGOLiN").info(f"Added info to {annotated} RGPs")

    def _write_graphs(self, output: Path, basename: str, graph_formats: list[str]):
        if "gexf" in graph_formats:
            graph_filename = output / f"{basename}.gexf"
            logging.getLogger("PPanGGOLiN").info(f"Writing RGP graph in GEXF format to {graph_filename}")
            nx.write_gexf(self.graph, graph_filename)

        if "graphml" in graph_formats:
            graph_filename = output / f"{basename}.graphml"
            logging.getLogger("PPanGGOLiN").info(f"Writing RGP graph in GraphML format to {graph_filename}")
            nx.write_graphml(self.graph, graph_filename)

    def _write_cluster_table(self, output: Path, basename: str, metric: str):
        outfile = output / f"{basename}.tsv"
        rows = []

        for node, attrs in self.graph.nodes(data=True):
            cluster = attrs.get(f"{metric}_cluster")
            if cluster is None:
                raise ValueError(
                    f"Node {node} does not have a '{metric}_cluster' attribute"
                )

            name = attrs["name"]

            if name.startswith("identical_rgps_"):
                identical_rgp_names = attrs["identical_rgp_names"]

                for child_name in identical_rgp_names.split(";"):
                    rows.append(
                        {
                            "RGPs": child_name,
                            "cluster": cluster,
                            "spot_id": self._rgp_to_spot.get(child_name, "No spot"),
                        }
                    )
            else:

                rows.append(
                    {
                        "RGPs": attrs.get("name", str(node)),
                        "cluster": cluster,
                        "spot_id": attrs.get("spot_id", "No spot"),
                    }
                )

        rows.sort(key=lambda row: (row["cluster"], row["RGPs"]))
        pd.DataFrame(rows, columns=["RGPs", "cluster", "spot_id"]).to_csv(
            outfile,
            sep="\t",
            index=False,
        )
        logging.getLogger("PPanGGOLiN").info(f"Writing RGP clusters in TSV format to {outfile}")

    def _write_outputs(
        self, output: Path, basename: str, graph_formats: list[str], metric: str
    ):
        self._write_graphs(output, basename, graph_formats)
        self._write_cluster_table(output, basename, metric)

    def run(self, options: RGPClusteringOptions):

        self._construct_regions(
            with_metadata=options.with_metadata,
            metadata_sources=options.metadata_sources or None,
            ignore_incomplete_rgp=options.ignore_incomplete_rgp,
        )

        self._compute_all_metrics(
            options.grr_cutoff, options.metric, disable_bar=options.disable_prog_bar
        )

        if options.unmerge_identical_rgps:
            self._add_edges_to_identical_rgps()

        self._louvain_clustering(options.metric)


        self._add_info_to_rgps()

        self._write_outputs(
            options.output, options.basename, options.graph_formats, options.metric
        )

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

def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provided by user
    """

    mk_outdir(args.output, args.force)

    if args.metadata_sep is not None:
        logging.getLogger("PPanGGOLiN").warning(
            "--metadata_sep is obsolete and ignored; metadata values are JSON-encoded in graphs."
        )

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
            metadata_sources=args.metadata_sources,
            ignore_incomplete_rgp=args.ignore_incomplete_rgp,
            disable_prog_bar=args.disable_prog_bar,
        )
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : Sub_parser for cluster_rgp command

    :return : Parser arguments for cluster_rgp command
    """
    parser = sub_parser.add_parser(
        "rgp_cluster", formatter_class=argparse.RawTextHelpFormatter
    )
    parser.description = "Cluster RGPs based on their gene families."
    parser.category = "Regions of Genomic Plasticity"
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
        default=["graphml"],
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
        default=None,
        help=argparse.SUPPRESS,
    )
