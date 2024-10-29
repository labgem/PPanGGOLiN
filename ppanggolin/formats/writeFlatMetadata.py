#!/usr/bin/env python3

# default libraries
import argparse
from collections import defaultdict
import logging
from typing import List
from pathlib import Path

# installed libraries
import pandas as pd

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir
from ppanggolin.formats.readBinaries import check_pangenome_info


def write_flat_metadata_files(
    pangenome: Pangenome,
    output: Path,
    pangenome_elements: List[str] = None,
    metadata_sources: List[str] = None,
    compress: bool = False,
    disable_bar: bool = False,
) -> None:
    """
    Main function to write flat metadata files from a pangenome.
    :todo: Split the function in subfunction

    :param pangenome: Pangenome object
    :param output: Path to output directory
    :param pangenome_elements: List of pangenome elements to include metadata for
    :param metadata_sources: List of metadata sources to use; if None, all sources are used
    :param compress: Compress the output files in .gz format
    :param disable_bar: Disable the progress bar
    """

    if not pangenome.has_metadata():
        logging.getLogger("PPanGGOLiN").warning(
            "No metadata is assigned to any pangenome element. Writing metadata is not possible."
        )
        return

    if pangenome_elements is None:
        pangenome_elements = []

    if all(
        pangenome.status["metadata"][element] == "No" for element in pangenome_elements
    ):
        logging.getLogger("PPanGGOLiN").warning(
            f"No metadata is assigned to any of the requested pangenome elements: {', '.join(pangenome_elements)}."
        )
        return

    if metadata_sources is None:
        # Use all available sources if none are specified
        element_to_sources = {
            element: sources
            for element, sources in pangenome.status["metasources"].items()
            if sources
        }
    else:
        # Use only the specified sources, if they match available sources
        element_to_sources = {}
        for element in pangenome_elements:
            existing_sources = pangenome.status["metasources"][element]
            matching_sources = set(existing_sources) & set(metadata_sources)
            if matching_sources:
                element_to_sources[element] = matching_sources

        if not element_to_sources:
            logging.getLogger("PPanGGOLiN").warning(
                f"None of the specified metadata sources ({metadata_sources}) match the requested pangenome elements: {pangenome_elements}."
            )
            return

    logging.getLogger("PPanGGOLiN").info(
        f"Writing metadata for {', '.join(element_to_sources.keys())} from {len([s for sources in element_to_sources.values() for s in sources])} sources."
    )

    need_dict = {
        "need_annotations": True,
        "need_families": "families" in element_to_sources,
        "need_rgp": "RGPs" in element_to_sources,
        "need_spots": "spots" in element_to_sources,
        "need_modules": "modules" in element_to_sources,
        "need_metadata": True,
        "sources": metadata_sources,
    }

    check_pangenome_info(pangenome, disable_bar=disable_bar, **need_dict)

    element_type_to_attribute = {
        "families": "gene_families",
        "genomes": "organisms",
        "RGPs": "regions",
        "genes": "genes",
        "modules": "modules",
        "spots": "spots",
        "contigs": "contigs",
    }

    for element_type, sources in element_to_sources.items():
        first_columns = [element_type]
        source_to_metadata = defaultdict(list)

        for element in getattr(pangenome, element_type_to_attribute[element_type]):
            for source in set(element.sources) & set(sources):
                for metadata in element.get_metadata_by_source(source).values():
                    metadata_dict = metadata.to_dict()
                    if element_type == "spots":
                        metadata_dict[element_type] = f"spot_{element.ID}"
                    elif element_type == "modules":
                        metadata_dict[element_type] = f"module_{element.ID}"
                    elif element_type in ["genes"]:
                        metadata_dict[element_type] = element.ID
                    elif element_type in ["families", "genomes", "RGPs", "contigs"]:
                        metadata_dict[element_type] = element.name

                    # Add genome name for genome-specific elements if genome is not already a metadata
                    if (
                        element_type in ["genes", "contigs", "RGPs"]
                        and "genomes" not in metadata_dict
                    ):
                        metadata_dict["genomes"] = element.organism.name
                        first_columns = ["genomes", element_type]

                    source_to_metadata[source].append(metadata_dict)

        for source, metadata_list in source_to_metadata.items():
            df_source = pd.DataFrame(metadata_list)
            columns_order = first_columns + [
                col for col in df_source.columns if col not in first_columns
            ]
            df_source = df_source.reindex(columns_order, axis=1)
            tsv_file_name = f"{element_type}_metadata_from_{source}.tsv"
            if compress:
                tsv_file_name += ".gz"
            tsv_path = output / tsv_file_name

            logging.getLogger("PPanGGOLiN").info(
                f"Writing {element_type} metadata from source '{source}' to '{tsv_path}'"
            )
            df_source.to_csv(tsv_path, sep="\t", index=False)

    logging.getLogger("PPanGGOLiN").info("Finished writing metadata in TSV format.")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    write_flat_metadata_files(
        pangenome,
        metadata_sources=args.metadata_sources,
        pangenome_elements=args.pangenome_elements,
        output=args.output,
        compress=args.compress,
        disable_bar=args.disable_prog_bar,
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "write_metadata", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_flat(parser)
    return parser


def parser_flat(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of write command

    :param parser: parser for align argument
    """
    pangenome_elements = {
        "families",
        "genomes",
        "contigs",
        "genes",
        "RGPs",
        "spots",
        "modules",
    }
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
        "--compress",
        required=False,
        action="store_true",
        help="Compress the files in .gz",
    )

    optional.add_argument(
        "-e",
        "--pangenome_elements",
        required=False,
        nargs="+",
        choices=pangenome_elements,
        default=pangenome_elements,
        help="Specify pangenome elements for which to write metadata. "
        "default is all element with metadata. ",
    )

    optional.add_argument(
        "-s",
        "--metadata_sources",
        default=None,
        nargs="+",
        help="Which source of metadata should be written. "
        "By default all metadata sources are included.",
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
