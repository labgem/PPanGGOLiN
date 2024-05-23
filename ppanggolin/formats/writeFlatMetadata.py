#!/usr/bin/env python3
# coding:utf-8

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
from ppanggolin.utils import write_compressed_or_not, mk_outdir, extract_contig_window, parse_input_paths_file
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.formats.write_proksee import write_proksee_organism
from ppanggolin.formats.writeSequences import read_genome_file, write_spaced_fasta
from ppanggolin.formats.writeFlatGenomes import get_organism_list 


def write_flat_metadata_files(pangenome: Pangenome, output: Path, 
                              pangenome_elements: List[str] = None, metadata_sources: List[str] = None, 
                               compress: bool = False, disable_bar: bool = False):
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

    if not pangenome.has_metadata():
        logging.getLogger("PPangGGOLiN").warning("No metadata is assigned to any pangenome element. Therefore it is not possible to write them.")


    elif all(( pangenome.status['metadata'][element] == "No" for element in pangenome_elements)):
            logging.getLogger("PPangGGOLiN").warning(f"No metadata is assigned to any of the requested pangenome elements: {','.join(pangenome_elements)}.")
            return

    if metadata_sources is None: 
        # source has not been specified and therefore all sources are considered
        element_to_sources = {element:sources  for element, sources in pangenome.status['metasources'].items() if sources}

    else:
        # sources has been specified. Checking that they match actual sources
        element_to_sources = {}
        for element in pangenome_elements:
            existing_sources = pangenome.status['metasources'][element]
            if set(existing_sources) & set(metadata_sources):
                element_to_sources[element] = set(existing_sources) & set(metadata_sources)

        if not element_to_sources:
            logging.getLogger("PPangGGOLiN").warning(f"None of given metadata sources ({metadata_sources}) "
                                                        f"is source of any of the requested pangenome elements: {pangenome_elements}.")
            return
        
    logging.getLogger("PPangGGOLiN").info(f"Writing metadata for {', '.join(element_to_sources)} from {len([s for sources in element_to_sources.values() for s in sources])} sources.")
    
    need_dict = {"need_annotations": True,
                 "need_families": "families" in element_to_sources,
                 "need_rgp": "RGPs" in element_to_sources,
                 "need_spots": "spots" in element_to_sources,
                 "need_modules": "modules" in element_to_sources, 
                 "need_metadata": True,
                 "sources": metadata_sources
                 }
    
    check_pangenome_info(pangenome, disable_bar=disable_bar, **need_dict,)

    element_type_to_attribute =  {"families":"gene_families", "genomes":"organisms", 
                           "RGPs":"regions", "genes":"genes", "modules":"modules",
                           "spots":"spots", "contigs":"contigs"}
    
    for element_type, sources in element_to_sources.items():
        first_columns = [element_type]
        source_to_metadata = defaultdict(list)

        for element in getattr(pangenome, element_type_to_attribute[element_type]):

            for source in set(element.sources) & set(sources):
            
                for metadata in element.get_metadata_by_source(source):
                    metadata_dict =  metadata.to_dict()
                    if element_type == "spots":
                        metadata_dict[element_type] = f"spot_{element.ID}"
                    
                    elif element_type == "modules":
                        metadata_dict[element_type] = f"module_{element.ID}"

                    elif element_type in ["genes", "contigs"]:
                        metadata_dict[element_type] = element.ID

                    elif element_type in ["families", 'genomes', "RGPs"]:
                        metadata_dict[element_type] = element.name

                    # add genome name for genome specific element if genomes is not a metadata
                    if element_type in ["genes", "contigs", "RGPs"] and "genomes" not in metadata_dict:
                        metadata_dict["genomes"] = element.organism.name
                        first_columns = ["genomes"] + first_columns

                    source_to_metadata[source].append(metadata_dict)


        for source, metadata_list in source_to_metadata.items():

            df_source = pd.DataFrame(metadata_list)
            columns_order = first_columns + [col for col in df_source.columns if col not in first_columns]  
            df_source = df_source.reindex(columns_order, axis=1)
            tsv_file_name = f"{element_type}_metadata_from_{source}.tsv" 
            if compress:
                tsv_file_name = f"{tsv_file_name}.gz"

            tsv_path = output / tsv_file_name

            logging.getLogger("PPangGGOLiN").info(f"Writing {element_type} metadata from source {source} in '{tsv_path}'")
            df_source.to_csv(tsv_path, sep="\t", index=False)

    logging.getLogger("PPangGGOLiN").info("Done writing metadata in TSV format")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    write_flat_metadata_files(pangenome, metadata_sources=args.metadata_sources,
                              pangenome_elements=args.pangenome_elements,
                               output=args.output, compress=args.compress,  disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("write_metadata", formatter_class=argparse.RawTextHelpFormatter)
    parser_flat(parser)
    return parser


def parser_flat(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of write command

    :param parser: parser for align argument
    """
    pangenome_elements = {"families", "genomes", "contigs", "genes", "RGPs", "spots", "modules"}
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=Path,
                          help="Output directory where the file(s) will be written")
    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument("--compress", required=False, action="store_true",
                          help="Compress the files in .gz")

    optional.add_argument("-e", "--pangenome_elements",
                          required=False,
                          nargs="+",
                          choices=pangenome_elements,
                          default=pangenome_elements,
                          help="Specify pangenome elements for which to write metadata. "
                               "default is all element with metadata. ")

    optional.add_argument("-s", "--metadata_sources",
                          default=None,
                          nargs="+",
                          help="Which source of metadata should be written. "
                               "By default all metadata sources are included.")


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
