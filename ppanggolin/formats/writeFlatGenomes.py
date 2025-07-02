#!/usr/bin/env python3

# default libraries
import argparse
from itertools import combinations
from collections import defaultdict
import logging
from typing import List, Dict, Set, Tuple
from pathlib import Path
import random
from tqdm import tqdm
import time
from concurrent.futures import ThreadPoolExecutor

# installed libraries
import networkx as nx
from plotly.express.colors import qualitative
import pandas as pd

# local libraries
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.genome import Organism, Gene, RNA
from ppanggolin.region import Region, Module
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import (
    write_compressed_or_not,
    mk_outdir,
    extract_contig_window,
    parse_input_paths_file,
)
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.formats.write_proksee import write_proksee_organism
from ppanggolin.formats.writeSequences import read_genome_file, write_spaced_fasta


def count_neighbors_partitions(gene_family: GeneFamily):
    """
    Count partition of neighbors families.

    :param gene_family: Gene family for which we count neighbors
    """

    nb_pers = 0
    nb_shell = 0
    nb_cloud = 0

    for neighbor in gene_family.neighbors:
        if neighbor.named_partition == "persistent":
            nb_pers += 1
        elif neighbor.named_partition == "shell":
            nb_shell += 1
        else:
            nb_cloud += 1

    return nb_pers, nb_shell, nb_cloud


def write_tsv_genome_file(
    organism: Organism,
    output: Path,
    compress: bool = False,
    metadata_sep: str = "|",
    need_regions: bool = False,
    need_spots: bool = False,
    need_modules: bool = False,
):
    """
    Write the table of genes with pangenome annotation for one organism in tsv

    :param organism: An organism
    :param output: Path to output directory
    :param compress: Compress the file in .gz
    :param need_regions: Write information about regions
    :param need_spots: Write information about spots
    :param need_modules: Write information about modules

    """

    rows = []
    for gene in organism.genes:
        nb_pers, nb_shell, nb_cloud = count_neighbors_partitions(gene.family)

        gene_info = {}

        gene_info["gene"] = (
            gene.ID if gene.local_identifier == "" else gene.local_identifier
        )
        gene_info["contig"] = gene.contig.name
        gene_info["start"] = gene.start
        gene_info["stop"] = gene.stop
        gene_info["strand"] = gene.strand
        gene_info["family"] = gene.family.name
        gene_info["nb_copy_in_genome"] = len(
            list(gene.family.get_genes_per_org(organism))
        )
        gene_info["partition"] = gene.family.named_partition
        gene_info["persistent_neighbors"] = nb_pers
        gene_info["shell_neighbors"] = nb_shell
        gene_info["cloud_neighbors"] = nb_cloud

        if need_regions:
            gene_info["RGP"] = str(gene.RGP) if gene.RGP is not None else None
        if need_spots:
            gene_info["Spot"] = str(gene.spot) if gene.spot is not None else None
        if need_modules:
            gene_info["Module"] = (
                str(gene.family.module) if gene.family.has_module else None
            )

        # Add metadata
        gene_metadata = {
            f"gene_{key}": value
            for key, value in gene.formatted_metadata_dict_to_string(
                metadata_sep
            ).items()
        }
        gene_info.update(gene_metadata)

        family_metadata = {
            f"family_{key}": value
            for key, value in gene.family.formatted_metadata_dict_to_string(
                metadata_sep
            ).items()
        }
        gene_info.update(family_metadata)

        if need_regions and gene.RGP:
            rgp_metadata = {
                f"rgp_{key}": value
                for key, value in gene.RGP.formatted_metadata_dict_to_string(
                    metadata_sep
                ).items()
            }
            gene_info.update(rgp_metadata)

        rows.append(gene_info)

    pd.DataFrame(rows).to_csv(
        output / f"{organism.name}.tsv{'.gz' if compress else ''}",
        sep="\t",
        index=False,
    )

    logging.getLogger("PPangGGOLiN").debug(
        f"Done writing the table with pangenome annotation for {organism.name}"
    )


def manage_module_colors(
    modules: Set[Module], window_size: int = 100
) -> Dict[Module, str]:
    """
    Manages colors for a list of modules based on gene positions and a specified window size.

    :param modules: A list of module objects for which you want to determine colors.
    :param window_size: Minimum number of genes between two modules to color them with the same color.
                        A higher value results in more module colors.
    :return: A dictionary that maps each module to its assigned color.
    """
    random.seed(1)

    color_mod_graph = nx.Graph()
    color_mod_graph.add_nodes_from(module for module in modules)

    contig_to_mod_genes = defaultdict(set)
    gene_to_module = {}

    for module in modules:
        for fam in module.families:
            for gene in fam.genes:
                contig_to_mod_genes[gene.contig].add(gene)
                gene_to_module[gene] = module

    for contig, mod_genes in contig_to_mod_genes.items():
        gene_positions = (gene.position for gene in mod_genes)
        contig_windows = extract_contig_window(
            contig.number_of_genes,
            gene_positions,
            window_size=window_size,
            is_circular=contig.is_circular,
        )
        contig_windows = list(contig_windows)

        for start, end in contig_windows:
            module_in_window = {
                gene_to_module[gene]
                for gene in mod_genes
                if start <= gene.position <= end
            }

            # Add edges between closely located modules
            module_edges = [
                (mod_a, mod_b) for mod_a, mod_b in combinations(module_in_window, 2)
            ]
            color_mod_graph.add_edges_from(module_edges)

    module_to_group = nx.coloring.greedy_color(color_mod_graph)

    # Attempt to have always the same color associated with the same module...
    module_to_color_int = {}
    group_with_color = []
    for module in sorted(modules, key=lambda x: x.ID):
        group = module_to_group[module]
        if group not in group_with_color:
            group_with_color.append(group)
        module_to_color_int[module] = group_with_color.index(group)

    # If you want to export the graph to see the coloring:
    # nx.set_node_attributes(color_mod_graph, module_to_color_int, name="color")
    # nx.readwrite.graphml.write_graphml(color_mod_graph, f"module_graph_window_size{window_size}.graphml")

    nb_colors = len(set(module_to_color_int.values()))
    logging.getLogger().debug(
        f"We have found that {nb_colors} colors were necessary to color Modules."
    )
    colors = palette(nb_colors)
    module_to_color = {mod: colors[col_i] for mod, col_i in module_to_color_int.items()}

    return module_to_color


def palette(nb_colors: int) -> List[str]:
    """
    Generates a palette of colors for visual representation.

    :param nb_colors: The number of colors needed in the palette.

    :return: A list of color codes in hexadecimal format.
    """

    # Combine two sets of predefined colors for variety
    colors = qualitative.Safe + qualitative.Vivid

    if len(colors) < nb_colors:
        # Generate random colors if not enough predefined colors are available
        random_colors = [
            "#" + "".join([random.choice("0123456789ABCDEF") for _ in range(6)])
            for _ in range(nb_colors - len(colors))
        ]
        colors += random_colors
    else:
        colors = colors[:nb_colors]

    return colors


def encode_attribute_val(product: str) -> str:
    """
    Encode special characters forbidden in column 9 of the GFF3 format.

    :param product: The input string to encode.
    :return: The encoded string with special characters replaced.

    Reference:
    - GFF3 format requirement: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    - Code source taken from Bakta: https://github.com/oschwengers/bakta
    """
    product = str(product)
    product = product.replace("%", "%25")
    product = product.replace(";", "%3B")
    product = product.replace("=", "%3D")
    product = product.replace("&", "%26")
    product = product.replace(",", "%2C")
    return product


def encode_attributes(attributes: List[Tuple]) -> str:
    """
    Encode a list of attributes in GFF3 format.

    :param attributes: A list of attribute key-value pairs represented as tuples.
    :return: The encoded attributes as a semicolon-separated string.
    """
    return ";".join(
        [
            f"{encode_attribute_val(k)}={encode_attribute_val(v)}"
            for k, v in attributes
            if str(v) != "" and v is not None
        ]
    )


def write_gff_file(
    organism: Organism,
    outdir: Path,
    annotation_sources: Dict[str, str],
    genome_sequences: Dict[str, str],
    metadata_sep: str = "|",
    compress: bool = False,
):
    """
    Write the GFF file of the provided organism.

    :param organism: Organism object for which the GFF file is being written.
    :param outdir: Path to the output directory where the GFF file will be written.
    :param metadata_sep: The separator used to join multiple metadata values
    :param compress: If True, compress the output GFF file using .gz format.
    :param annotation_sources: A dictionary that maps types of features to their source information.
    :param genome_sequences: A dictionary mapping contig names to their DNA sequences (default: None).
    """

    # sort contig by their name
    sorted_contigs = sorted(organism.contigs, key=lambda x: x.name)

    organism_metadata = [
        (f"genome_{key}", value)
        for key, value in organism.formatted_metadata_dict_to_string(
            metadata_sep
        ).items()
    ]

    with write_compressed_or_not(outdir / f"{organism.name}.gff", compress) as outfile:
        # write gff header
        outfile.write("##gff-version 3\n")
        for contig in sorted_contigs:
            if contig.length is None:
                raise AttributeError(f"Contig {contig.name} has no length defined.")

            outfile.write(f"##sequence-region {contig.name} 1 {contig.length}\n")

        for contig in sorted_contigs:
            contig_metadata = [
                (f"contig_{key}", value)
                for key, value in contig.formatted_metadata_dict_to_string(
                    metadata_sep
                ).items()
            ]
            attributes = (
                [
                    ("ID", contig.name),
                    ("Is_circular", "true" if contig.is_circular else "false"),
                ]
                + organism_metadata
                + contig_metadata
            )
            attributes_str = encode_attributes(attributes)

            contig_line = [
                contig.name,
                ".",
                "region",
                "1",
                contig.length,
                ".",
                "+",
                ".",
                attributes_str,
            ]
            contig_line_str = "\t".join(map(str, contig_line))
            outfile.write(contig_line_str + "\n")

            contig_elements = sorted(
                list(contig.regions) + list(contig.genes) + list(contig.RNAs),
                key=lambda x: x.start,
            )

            for feature in contig_elements:
                phase = "."

                if isinstance(feature, (Gene, RNA)):
                    feat_type = feature.type

                    strand = feature.strand

                    source = annotation_sources.get(feat_type, "external")

                    # before the CDS or RNA line a gene line is created. with the following id
                    parent_gene_id = f"gene-{feature.ID}"

                    attributes = [
                        ("ID", feature.ID),
                        ("Name", feature.name),
                        ("Parent", parent_gene_id),
                        ("product", feature.product),
                    ]

                    score = "."

                    if isinstance(feature, Gene):
                        rgp = feature.RGP.name if feature.RGP else ""
                        phase = "0"

                        attributes += [
                            ("family", feature.family.name),
                            ("partition", feature.family.named_partition),
                            ("rgp", rgp),
                            (
                                "module",
                                feature.family.module,
                            ),  # family.module can be None...
                        ]

                        # adding attributes
                        gene_metadata = [
                            (f"gene_{key}", value)
                            for key, value in feature.formatted_metadata_dict_to_string(
                                metadata_sep
                            ).items()
                        ]
                        family_metadata = [
                            (f"family_{key}", value)
                            for key, value in feature.family.formatted_metadata_dict_to_string(
                                metadata_sep
                            ).items()
                        ]

                        attributes += gene_metadata
                        attributes += family_metadata

                    # add an extra line of type gene
                    stop = feature.stop
                    if feature.overlaps_contig_edge:
                        stop = contig.length + feature.stop

                    gene_line = [
                        contig.name,
                        source,
                        "gene",
                        feature.start,
                        stop,
                        ".",
                        strand,
                        ".",
                        f"ID={encode_attribute_val(parent_gene_id)}",
                    ]
                    line_str = "\t".join(map(str, gene_line))
                    outfile.write(line_str + "\n")

                elif isinstance(feature, Region):
                    feat_type = "region"
                    source = "ppanggolin"
                    strand = "."
                    score = "."

                    rgp_metadata = [
                        (f"rgp_{key}", value)
                        for key, value in feature.formatted_metadata_dict_to_string(
                            metadata_sep
                        ).items()
                    ]

                    attributes = [
                        ("Name", feature.name),
                        (
                            "spot",
                            feature.spot.ID if feature.spot is not None else "No_spot",
                        ),
                        ("Note", "Region of Genomic Plasticity (RGP)"),
                    ]
                    attributes += rgp_metadata

                else:
                    raise TypeError(
                        f"The feature to write in gff file does not have an expected types. {type(feature)}"
                    )

                attributes_str = encode_attributes(attributes)

                coordinates = feature.coordinates
                if feature.overlaps_contig_edge:
                    coordinates = convert_overlapping_coordinates_for_gff(
                        feature.coordinates, len(contig)
                    )

                for start, stop in coordinates:
                    line = [
                        contig.name,
                        source,  # Source
                        feat_type,
                        start,
                        stop,
                        score,
                        strand,
                        phase,
                        attributes_str,
                    ]

                    line_str = "\t".join(map(str, line))
                    outfile.write(line_str + "\n")

        if genome_sequences:
            logging.getLogger("PPanGGOLiN").debug(
                "Writing fasta section of gff file..."
            )
            outfile.write("##FASTA\n")
            for contig in sorted_contigs:
                outfile.write(f">{contig.name}\n")

                outfile.write(
                    write_spaced_fasta(genome_sequences[contig.name], space=60)
                )


def convert_overlapping_coordinates_for_gff(
    coordinates: List[Tuple[int, int]], contig_length: int
):
    """
    Converts overlapping gene coordinates in GFF format for circular contigs.

    :param coordinates: List of tuples representing gene coordinates.
    :param contig_length: Length of the circular contig.
    """

    start, stop = coordinates[0]
    new_coordinates = [(start, stop)]
    # convert all coordinates that are at the beginning
    # of the contig to the extent of the contig
    for start_n, stop_n in coordinates[1:]:
        if start_n < start:  # we are on the beginning of the contig
            new_start = contig_length + start_n
            new_stop = contig_length + stop_n
            new_coordinates.append((new_start, new_stop))
            start, stop = new_start, new_stop

        else:
            start, stop = start_n, stop_n

            new_coordinates.append((start, stop))

    # merge continuous coordinates
    merged_coordinates = []
    start, stop = new_coordinates[0]

    for start_n, stop_n in new_coordinates[1:]:
        if stop + 1 == start_n:
            stop = stop_n
        else:
            merged_coordinates.append((start, stop))
            start, stop = start_n, stop_n

    merged_coordinates.append((start, stop))

    return merged_coordinates


def get_organism_list(organisms_filt: str, pangenome: Pangenome) -> Set[Organism]:
    """
    Get a list of organisms to include in the output.

    :param organisms_filt: Filter for selecting organisms. It can be a file path with one organism name per line
                          or a comma-separated list of organism names.
    :param pangenome: The pangenome from which organisms will be selected.
    :return: A set of selected Organism objects.
    """
    if organisms_filt == "all":
        logging.getLogger("PPanGGOLiN").info(
            "Writing output for all genomes of the pangenome."
        )
        organisms_list = set(pangenome.organisms)

    else:
        if Path(organisms_filt).is_file():
            logging.getLogger("PPanGGOLiN").debug(
                "Parsing the list of genomes from a file "
                "to determine which genomes should be included in the output."
            )
            with open(organisms_filt) as fl:
                org_names = [
                    line.strip() for line in fl if line and not line.startswith("#")
                ]
        else:
            org_names = [
                name.strip() for name in organisms_filt.split(",") if name.strip()
            ]

        organisms_list = set()
        org_not_in_pangenome = set()
        for org_name in org_names:
            try:
                org = pangenome.get_organism(org_name)
                organisms_list.add(org)
            except KeyError:
                org_not_in_pangenome.add(org_name)
        if org_not_in_pangenome:
            raise KeyError(
                f"{len(org_not_in_pangenome)} organism(s) specified with '--genomes' parameter were "
                f"not found in the pangenome: {', '.join(org_not_in_pangenome)}"
            )

        logging.getLogger("PPanGGOLiN").info(
            f"Writing output for {len(organisms_list)}/{pangenome.number_of_organisms} genomes of the pangenome."
        )

    return organisms_list


def mp_write_genomes_file(
    organism: Organism,
    output: Path,
    organisms_file: Path = None,
    proksee: bool = False,
    gff: bool = False,
    table: bool = False,
    **kwargs,
) -> str:
    """Wrapper for the write_genomes_file function that allows it to be used in multiprocessing.

    :param organism: Specify the organism to be written
    :param output: Specify the path to the output directory
    :param organisms_file: Read the genome sequences from a file
    :param proksee: Write a proksee file for the organism
    :param gff:  Write the gff file for the organism
    :param table: Write the organism file for the organism
    :param kwargs: Pass any number of keyword arguments to the function

    :return: The organism name
    """

    genome_sequences = None
    if organisms_file and (gff or proksee):
        genome_sequences = read_genome_file(organisms_file, organism)

    if proksee:

        proksee_outdir = output / "proksee"
        mk_outdir(proksee_outdir, force=True, exist_ok=True)
        output_file = proksee_outdir / f"{organism.name}.json"

        # Write ProkSee data for the organism
        write_proksee_organism(
            organism,
            output_file,
            features=["all"],
            genome_sequences=genome_sequences,
            **{
                arg: kwargs[arg]
                for arg in kwargs.keys()
                & {"module_to_colors", "compress", "multigenics"}
            },
        )

    if gff:
        gff_outdir = output / "gff"
        mk_outdir(gff_outdir, force=True, exist_ok=True)

        write_gff_file(
            organism,
            outdir=gff_outdir,
            genome_sequences=genome_sequences,
            **{
                arg: kwargs[arg]
                for arg in kwargs.keys()
                & {"compress", "annotation_sources", "metadata_sep"}
            },
        )

    if table:
        table_outdir = output / "table"
        mk_outdir(table_outdir, force=True, exist_ok=True)

        write_tsv_genome_file(
            organism=organism,
            output=table_outdir,
            **{
                arg: kwargs[arg]
                for arg in kwargs.keys()
                & {
                    "need_regions",
                    "need_modules",
                    "need_spots",
                    "compress",
                    "metadata_sep",
                }
            },
        )

    return organism.name


def write_flat_genome_files(
    pangenome: Pangenome,
    output: Path,
    table: bool = False,
    gff: bool = False,
    proksee: bool = False,
    compress: bool = False,
    fasta: Path = None,
    anno: Path = None,
    organisms_filt: str = "all",
    add_metadata: bool = False,
    metadata_sep: str = "|",
    metadata_sources: List[str] = None,
    cpu: int = 1,
    disable_bar: bool = False,
):
    """
    Main function to write flat files from pangenome

    :param pangenome: Pangenome object
    :param output: Path to output directory
    :param cpu: Number of available core
    :param table: write table with pangenome annotation for each genome
    :param gff: write a gff file with pangenome annotation for each organism
    :param proksee: write a proksee file with pangenome annotation for each organisms
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    :param fasta: File containing the list FASTA files for each organism
    :param anno: File containing the list of GBFF/GFF files for each organism
    :param organisms_filt: String used to specify which organism to write. if all, all organisms are written.
    :param add_metadata: Add metadata to GFF files
    :param metadata_sep: The separator used to join multiple metadata values
    :param metadata_sources: Sources of the metadata to use and write in the outputs. None means all sources are used.
    """

    if not any(x for x in [table, gff, proksee]):
        raise argparse.ArgumentError(
            argument=None, message="You did not indicate what file you wanted to write."
        )

    need_dict = {
        "need_annotations": True,
        "need_families": True,
        "need_partitions": True,
        "need_rgp": True if pangenome.status["predictedRGP"] != "No" else False,
        "need_spots": True if pangenome.status["spots"] != "No" else False,
        "need_modules": True if pangenome.status["modules"] != "No" else False,
        "need_graph": True if table else False,
        "need_metadata": add_metadata,
        "sources": metadata_sources,
    }

    # Place here to raise an error if file doesn't found before to read pangenome
    organisms_file = fasta if fasta is not None else anno

    check_pangenome_info(pangenome, disable_bar=disable_bar, **need_dict)

    organisms_list = get_organism_list(organisms_filt, pangenome)
    if not organisms_list:
        raise ValueError(
            "No genomes are selected for output. Please check the '--genomes' parameter."
        )

    multigenics = None
    if need_dict["need_rgp"]:
        multigenics = pangenome.get_multigenics(
            pangenome.parameters["rgp"]["dup_margin"]
        )

    org_dict = (
        parse_input_paths_file(organisms_file)
        if organisms_file and (gff or proksee)
        else None
    )

    if proksee:
        # Generate a color mapping for modules
        module_to_colors = manage_module_colors(set(pangenome.modules))

    organism2args = defaultdict(
        lambda: {
            "output": output,
            "table": table,
            "gff": gff,
            "proksee": proksee,
            "compress": compress,
            "multigenics": multigenics,
        }
    )
    for organism in organisms_list:
        organism_args = {
            "genome_file": org_dict[organism.name]["path"] if org_dict else None,
            "metadata_sep": metadata_sep,
        }

        if proksee:
            organism_args["module_to_colors"] = {
                module: module_to_colors[module] for module in organism.modules
            }

        if gff:
            # prepare variable for gff output
            if pangenome.parameters["annotate"]["# read_annotations_from_file"]:
                organism_args["annotation_sources"] = {
                    "rRNA": "external",
                    "tRNA": "external",
                    "CDS": "external",
                }
            else:
                organism_args["annotation_sources"] = {}

        if table:
            # create _genePerOrg dict with get_org_dict methodbefore the multiprocessing to prevent putative errors.
            # As this is used in multiprocessing when computing nb_copy_in_genome.
            for family in pangenome.gene_families:
                family.get_org_dict()

            organism_args.update(
                {
                    "need_regions": need_dict["need_rgp"],
                    "need_modules": need_dict["need_modules"],
                    "need_spots": need_dict["need_spots"],
                }
            )
        organism2args[organism].update(organism_args)

    start_writing = time.time()
    with ThreadPoolExecutor(max_workers=cpu) as executor:
        with tqdm(
            total=(len(organisms_list)), unit="genome", disable=disable_bar
        ) as progress:
            futures = []
            for organism, kwargs in organism2args.items():
                logging.getLogger("PPanGGOLiN").debug(
                    f"Writing genome annotations for {organism.name}"
                )
                future = executor.submit(mp_write_genomes_file, organism, **kwargs)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                org_name = future.result()
                logging.getLogger("PPanGGOLiN").debug(
                    f"Done writing the GFF file with pangenome annotation for {org_name}."
                )

    writing_time = time.time() - start_writing
    logging.getLogger("PPanGGOLiN").debug(
        f"writing_time for {pangenome.number_of_organisms} genomes: {writing_time} seconds"
    )


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    write_flat_genome_files(
        pangenome,
        args.output,
        table=args.table,
        gff=args.gff,
        proksee=args.proksee,
        compress=args.compress,
        fasta=args.fasta,
        anno=args.anno,
        organisms_filt=args.genomes,
        add_metadata=args.add_metadata,
        metadata_sep=args.metadata_sep,
        metadata_sources=args.metadata_sources,
        cpu=args.cpu,
        disable_bar=args.disable_prog_bar,
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "write_genomes", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_flat(parser)
    return parser


def parser_flat(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of write command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome .h5 file"
    )
    required.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output directory where the file(s) will be written",
    )
    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument(
        "--table",
        required=False,
        action="store_true",
        help="Generate a tsv file for each genome with pangenome annotations.",
    )

    optional.add_argument(
        "--gff",
        required=False,
        action="store_true",
        help="Generate a gff file for each genome containing pangenome annotations.",
    )

    optional.add_argument(
        "--proksee",
        required=False,
        action="store_true",
        help="Generate JSON map files for PROKSEE for each genome containing pangenome annotations "
        "to be used in proksee.",
    )

    optional.add_argument(
        "--compress",
        required=False,
        action="store_true",
        help="Compress the files in .gz",
    )

    optional.add_argument(
        "--genomes",
        required=False,
        default="all",
        help="Specify the genomes for which to generate output. "
        "You can provide a list of genome names either directly in the command line separated "
        "by commas, or by referencing a file containing the list of genome names, "
        "with one name per line.",
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

    optional.add_argument(
        "-c",
        "--cpu",
        required=False,
        default=1,
        type=int,
        help="Number of available cpus",
    )

    context = parser.add_argument_group(
        title="Contextually required arguments",
        description="With --proksee and --gff, the following arguments can be "
        "used to add sequence information to the output file:",
    )

    context.add_argument(
        "--fasta",
        required=False,
        type=Path,
        help="A tab-separated file listing the genome names, and the fasta filepath of its genomic "
        "sequence(s) (the fastas can be compressed with gzip). One line per genome.",
    )

    context.add_argument(
        "--anno",
        required=False,
        type=Path,
        help="A tab-separated file listing the genome names, and the gff/gbff filepath of its "
        "annotations (the files can be compressed with gzip). One line per genome. "
        "If this is provided, those annotations will be used.",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_flat(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
