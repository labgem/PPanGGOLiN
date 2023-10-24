#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from multiprocessing import get_context
from itertools import combinations
from collections import Counter, defaultdict
import logging
from typing import TextIO,List, Dict
from pathlib import Path
from typing import TextIO
from importlib.metadata import distribution
from statistics import median, mean, stdev
import csv
import random
from tqdm import tqdm
import time

import networkx as nx
from plotly.express.colors import qualitative


# local libraries
from ppanggolin.edge import Edge
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.genome import Organism, Gene, Contig, RNA
from ppanggolin.region import Region, Spot, Module
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import write_compressed_or_not, mk_outdir, restricted_float, extract_contig_window, parse_input_paths_file
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.formats.write_proksee import write_proksee_organism
from ppanggolin.formats.writeSequences import read_genome_file, write_spaced_fasta

def count_neighbors_partitions(gene_family:GeneFamily):
    """
    Count partition of neighbors families.

    :param gene_family: Gene family forwhich we count neighbors

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

def write_org_file(org: Organism, output: Path, compress: bool = False, 
                   add_regions:bool = False, add_spots:bool = False, add_modules:bool = False):
    """
    Write the table of pangenome for one organism in tsv format

    :param org: Projected organism
    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """

    with write_compressed_or_not(output / f"{org.name}.tsv", compress) as outfile:
        fieldnames = ["gene", "contig", "start", "stop", "strand", "family", "nb_copy_in_org",
                  "partition", "persistent_neighbors", "shell_neighbors", "cloud_neighbors"]
        
        if add_regions:
            fieldnames.append("RGPs")
        if add_spots:
            fieldnames.append("Spots")
        if add_modules:
            fieldnames.append("Modules")

        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for contig in org.contigs:
            for gene in contig.genes:
                
                nb_pers,nb_shell,nb_cloud = count_neighbors_partitions(gene.family)

                row = {"gene":gene.ID if gene.local_identifier == "" else gene.local_identifier,
                       "contig":contig.name,
                       "start":gene.start,
                       "stop":gene.stop,
                       "strand":gene.strand,
                       "family":gene.family.name,
                       "nb_copy_in_org":len(list(gene.family.get_genes_per_org(org))), 
                       "partition":gene.family.named_partition,
                       "persistent_neighbors":nb_pers,
                       "shell_neighbors":nb_shell,
                       "cloud_neighbors":nb_cloud}
                
                if add_regions:
                    row["RGPs"] = gene.RGP.name if gene.RGP is not None else None

                if add_spots:
                    spot = None
                    if gene.family.number_of_spots > 0:
                        spot = ','.join([str(spot.ID) for spot in gene.family.spots])
                    row['Spots'] = spot

                if add_modules:
                    modules = None
                    if gene.family.number_of_modules > 0:
                        modules = ','.join(["module_" + str(module.ID) for module in gene.family.modules])
                    row['Modules'] = modules

                writer.writerow(row)
        
    logging.getLogger().debug(f"Done writing the table with pangenome annotation for {org.name}")


def manage_module_colors(modules: List[Module], window_size:int=50) -> Dict[Module, str]:
    """
    Manages colors for a list of modules based on gene positions and a specified window size.

    :param modules: A list of module objects for which you want to determine colors.
    :param window_size: Minimum number of genes between two modules to color them with the same color. 
                        A higher value results in more module colors.
    :return: A dictionary that maps each module to its assigned color.
    """
    
    color_mod_graph = nx.Graph()
    color_mod_graph.add_nodes_from((module for module in modules))

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
            contig.number_of_genes, gene_positions, window_size=window_size, is_circular=contig.is_circular
        )
        contig_windows = list(contig_windows)

        for (start, end) in contig_windows:
            module_in_window = {gene_to_module[gene] for gene in mod_genes if start <= gene.position <= end}

            # Add edges between closely located modules
            module_edges = [(mod_a, mod_b) for mod_a, mod_b in combinations(module_in_window, 2)]
            color_mod_graph.add_edges_from(module_edges)

    module_to_color_int = nx.coloring.greedy_color(color_mod_graph)

    # If you want to export the graph to see the coloring:
    # nx.set_node_attributes(color_mod_graph, color_dict, name="color")
    # nx.readwrite.graphml.write_graphml(color_mod_graph, f"module_graph_window_size{window_size}.graphml")
    
    nb_colors = len(set(module_to_color_int.values()))
    logging.getLogger().debug(f"We have found that {nb_colors} colors were necessary to color Modules.")
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
    colors = qualitative.Vivid + qualitative.Safe
    
    if len(colors) < nb_colors:
        # Generate random colors if not enough predefined colors are available
        random.seed(1)
        random_colors = ["#" + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)]) for _ in range(nb_colors - len(colors))]
        colors += random_colors
    else:
        colors =  colors[:nb_colors]

    return colors


def write_gff_file(org: Organism, contig_to_rgp: Dict[Contig, Region], 
                   rgp_to_spotid: Dict[Region, str], outdir: str, compress: bool,
                   annotation_sources: Dict[str, str], genome_sequences:Dict[str,str]):
    """
    Write the GFF file of the provided organism.

    :param org: Organism object for which the GFF file is being written.
    :param contig_to_rgp: Dictionary mapping Contig objects to their corresponding Region objects.
    :param rgp_to_spotid: Dictionary mapping Region objects to their corresponding spot IDs.
    :param outdir: Path to the output directory where the GFF file will be written.
    :param compress: If True, compress the output GFF file using .gz format.
    :param annotation_sources: A dictionary that maps types of features to their source information.
    :param genome_sequences: A dictionary mapping contig names to their DNA sequences (default: None).
    """

    # sort contig by their name
    sorted_contigs = sorted(org.contigs, key= lambda x : x.name)

    with write_compressed_or_not(outdir /  F"{org.name}.gff", compress) as outfile:
        # write gff header
        outfile.write('##gff-version 3\n')
        for contig in sorted_contigs:
            if contig.length is None:
                raise AttributeError(f'Contig {contig.name} has no length defined.')

            outfile.write(f'##sequence-region {contig.name} 1 {contig.length}\n')

        for contig in sorted_contigs:
            contig_elements = sorted(contig_to_rgp[contig] + list(contig.genes) + list(contig.RNAs), key=lambda x: (x.start))

            for feature in contig_elements:
                
                if type(feature) in [Gene, RNA]:
                    feat_type = feature.type

                    strand = feature.strand
                    
                    source = annotation_sources.get(feat_type, "external")

                    # before the CDS or RNA line a gene line is created. with the following id
                    parent_gene_id=f"gene-{feature.ID}"

                    attributes = [("ID", feature.ID), 
                                  ("Name", feature.name),
                                  ('Parent', parent_gene_id),
                                  ("product", feature.product),
                                ]
                    
                    score = '.'
                    
                    if type(feature) == Gene:
                        rgp = feature.RGP.name if feature.RGP else ""
                        attributes += [
                            ("Family", feature.family.name),
                            ("Partition", feature.family.named_partition),
                            ('RGP', rgp),
                            ('Module', ','.join((f"module_{module.ID}" for module in feature.family.modules)) )
                        ]
                
                    # add an extra line of type gene
                    gene_line = [contig.name,
                            source, 
                            'gene',
                            feature.start,
                            feature.stop,
                            '.',
                            strand,
                            ".",
                            f'ID={parent_gene_id}'
                            ]
                    
                    line_str = '\t'.join(map(str, gene_line))
                    outfile.write(line_str + "\n")

                elif type(feature) == Region:
                    feat_type = "region"
                    source = "ppanggolin"
                    strand = "."
                    score = feature.score # TODO does RGP score make sens and do we want it in gff file?
                    attributes = [
                            ("Name", feature.name),
                            ("Spot", rgp_to_spotid.get(feature, "No_spot")),
                            ("Note", "Region of Genomic Plasticity (RGP)")
                    ]

                
                else:
                    raise TypeError(f'The feature to write in gff file does not have an expected types. {type(feature)}')


                attributes_str = ';'.join([f"{k}={v}" for k,v in attributes if v != "" and v is not None])

                line = [contig.name,
                        source, # Source
                        feat_type,
                        feature.start,
                        feature.stop,
                        score,
                        strand,
                        ".",
                        attributes_str,
                        ]

                line_str = '\t'.join(map(str, line))
                outfile.write(line_str + "\n")

        if genome_sequences:
            logging.getLogger("PPanGGOLiN").debug("Writing fasta section of gff file...")
            outfile.write(f"##FASTA\n")
            for contig in sorted_contigs:
                outfile.write(f">{contig.name}\n")

                outfile.write(write_spaced_fasta(genome_sequences[contig.name], space=60))


def write_flat_genome_files(pangenome: Pangenome, output: Path,
                            table: bool = False, gff: bool = False, proksee: bool = False, compress: bool = False,
                     disable_bar: bool = False, fasta=None, anno=None):
    """
    Main function to write flat files from pangenome

    :param pangenome: Pangenome object
    :param output: Path to output directory
    :param cpu: Number of available core
    :param soft_core: Soft core threshold to use
    :param dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated
    :param table: write table with pangenome annotation for each genome
    :param gff: write a gff file with pangenome annotation for each organisms
    :param proksee: write a proksee file with pangenome annotation for each organisms
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """

    if not any(x for x in [ table, gff, proksee]):
        raise Exception("You did not indicate what file you wanted to write.")

    needAnnotations = True
    needFamilies = True
    needPartitions = True
    needRegions = True if pangenome.status["predictedRGP"] == "inFile" else False
    needSpots = True if pangenome.status["spots"] == "inFile" else False
    needModules = True if pangenome.status["modules"] == "inFile" else False
    if table:
        need_graph = True
    else:
        need_graph = False

    # TODO: see if we export metadata in write_genomes command
    check_pangenome_info(pangenome, need_annotations=needAnnotations, need_families=needFamilies,need_graph=need_graph,
                         need_partitions=needPartitions, need_rgp=needRegions, need_spots=needSpots,
                         need_modules=needModules, need_metadata=None, metatype=None, sources=None,
                         disable_bar=disable_bar)
    
    organisms_file = fasta if fasta is not None else anno

    if organisms_file and (gff or proksee):
        org_dict = parse_input_paths_file(organisms_file)

    if proksee:
        org_to_modules = defaultdict(set)
        # Create a mapping of organisms and modules they have
        for mod in pangenome.modules:
            for org in mod.organisms:
                org_to_modules[org].add(mod)


    # Generate a color mapping for modules
    module_to_colors = manage_module_colors(set(pangenome.modules))


    if gff:
        # prepare variable for gff output
        if pangenome.parameters["annotate"]["# read_annotations_from_file"]:
            annotation_sources = {"rRNA": "external",
                                "tRNA": "external",
                                "CDS":"external"}
        else:
            annotation_sources = {}

        rgp_to_spot_id = {rgp:f"spot_{spot.ID}" for spot in pangenome.spots for rgp in spot.regions}

        contig_to_rgp = defaultdict(list)
        for rgp in pangenome.regions:
            contig_to_rgp[rgp.contig].append(rgp)
            rgp_to_spot_id = {rgp:f"spot_{spot.ID}" for spot in pangenome.spots for rgp in spot.regions}

    start_writing = time.time()
    
    # TODO try to multithread this part... ? 
    for organism in  tqdm(pangenome.organisms, total=pangenome.number_of_organisms, unit="organism", disable=disable_bar):

        logging.getLogger("PPanGGOLiN").debug(f"Writing genome annotations for {organism.name}")


        org_outdir = output / organism.name
        
        mk_outdir(org_outdir, force=True)

        genome_sequences = None
        if organisms_file and (gff or proksee):
            genome_sequences = read_genome_file(org_dict[organism.name]['path'], organism)

        if proksee:
            # Generate a color mapping for modules specific to the organism
            org_module_to_color = {org_mod: module_to_colors[org_mod] for org_mod in org_to_modules[organism]}

            output_file = org_outdir / f"{organism.name}.json"

            # TODO Add spot ID in attribute of RGP element in proksee

            # Write ProkSee data for the organism
            write_proksee_organism(organism, output_file, features=['all'], module_to_colors=org_module_to_color, rgps=pangenome.regions,
                                    genome_sequences=genome_sequences)
        
        
        if gff:
            write_gff_file(organism, contig_to_rgp, rgp_to_spot_id, org_outdir, compress, annotation_sources, genome_sequences)

        if table:
            write_org_file(org=organism, output=org_outdir, compress=compress,
                       add_regions=needRegions, add_modules=needModules, add_spots=needSpots)

    writing_time = time.time() - start_writing
    logging.getLogger("PPanGGOLiN").debug(f"writing_time for {pangenome.number_of_organisms} genomes: {writing_time} seconds")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    write_flat_genome_files(pangenome, args.output,
                    table=args.table, gff=args.gff, proksee=args.proksee,
                    compress=args.compress, disable_bar=args.disable_prog_bar, fasta=args.fasta, anno=args.anno)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("write_genomes", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_flat(parser)
    return parser


def parser_flat(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of write command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=Path,
                          help="Output directory where the file(s) will be written")
    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument("--table", required=False, action="store_true",
                          help="Generate a tsv file for each genome with pangenome annotations.")
    optional.add_argument("--gff", required=False, action="store_true",
                        help="Generate a gff file for each genome containing pangenome annotations.")
    
    optional.add_argument("--proksee", required=False, action="store_true",
                        help="Generate JSON map files for PROKSEE for each genome containing pangenome annotations to be used in proksee.")
    
    optional.add_argument("--compress", required=False, action="store_true", help="Compress the files in .gz")

    context = parser.add_argument_group(title="Contextually required arguments",
                                        description="With --proksee and -gff, the following arguments can be "
                                        "used to add sequence information to the output file:")
    
    context.add_argument('--fasta', required=False, type=Path,
                         help="A tab-separated file listing the organism names, and the fasta filepath of its genomic "
                              "sequence(s) (the fastas can be compressed with gzip). One line per organism.")
    
    context.add_argument('--anno', required=False, type=Path,
                         help="A tab-separated file listing the organism names, and the gff/gbff filepath of its "
                              "annotations (the files can be compressed with gzip). One line per organism. "
                              "If this is provided, those annotations will be used.")

if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_flat(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
