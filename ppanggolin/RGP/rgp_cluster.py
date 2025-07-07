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
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region, Spot, Module
from ppanggolin.formats import check_pangenome_info
from ppanggolin.utils import restricted_float, mk_outdir
from ppanggolin.geneFamily import GeneFamily


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
