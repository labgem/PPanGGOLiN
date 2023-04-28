#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
from pathlib import Path
import csv
import re

# installed libraries
from tqdm import tqdm
import pandas as pd

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.metadata import Metadata
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome


def check_pangenome_metadata(pangenome: Pangenome, source: str, metatype: str, force: bool = False,
                             disable_bar: bool = False):
    """ Check and load pangenome information before adding metadata

    :param pangenome: Pangenome object
    :param source: source of the metadata
    :param metatype: select to which pangenome element metadata will be added
    :param force: erase if a metadata for the provide source and metatype already exist
    :param disable_bar: Disable bar
    """
    need_dic = {'need_annotations': False,
                'need_families': False,
                'need_rgp': False,
                'need_spots': False,
                'need_modules': False}
    if metatype in ["genes", "genomes", "families"]:
        need_dic['need_annotations'] = True
    if metatype == "families":
        need_dic['need_families'] = True
    if metatype in ["RGPs", "spots"]:
        need_dic['need_rgp'] = True
    if metatype in ["RGPs", "spots", "modules"]:
        need_dic['need_spots'] = True
    if metatype in ["RGPs", "spots", "modules"]:
        need_dic['need_modules'] = True

    if pangenome.status["metadata"][metatype] == "inFile" and source in pangenome.status["metasources"][metatype]:
        if force:
            erase_pangenome(pangenome, metadata=True, source=source, metatype=metatype)
        else:
            raise Exception(
                f"An metadata corresponding to the source : '{source}' already exist in pangenome organims."
                "Add the option --force to erase")
    check_pangenome_info(pangenome, disable_bar=disable_bar, **need_dic)


def check_metadata_format(metadata: Path, metatype: str) -> pd.DataFrame:
    """Check if the TSV with metadata respect the input format

    :param metadata: Path to the TSV file with metadata
    :param metatype: Indicate which pangenome element metadata will be added

    :return: Dataframe with metadata loaded
    """
    assert metatype in ["families", "genomes", "genes", "RGPs", "spots", "modules"]

    colname_check = re.compile('^[a-zA-Z_][a-zA-Z0-9_]*$')
    metadata_df = pd.read_csv(metadata, sep="\t", header=0, quoting=csv.QUOTE_NONE,
                              dtype={metatype: str})
    metadata_df.replace(to_replace='-', value=pd.NA, inplace=True)
    if metatype not in metadata_df.columns or metadata_df.shape[1] < 2:
        raise KeyError(f"You should at least provide in columns names : {metatype} and one another value. "
                       "Look at documentation for more information")
    for column in metadata_df.columns:
        if not colname_check.match(column):
            raise ValueError(f"column name is not a valid identifier: {column}; "
                             f"it does not match the pattern {colname_check.pattern}")
        if column != metatype and metadata_df.dtypes[column] == object:
            pd.to_numeric(metadata_df[column], downcast='integer', errors='ignore')

    return metadata_df


def assign_metadata(metadata_df: pd.DataFrame, pangenome: Pangenome, source: str, metatype: str,
                    omit: bool = False, disable_bar: bool = False):
    """ Add to pangenome element a metadata

    :param metadata_df: Dataframe with for each family a metadata
    :param pangenome: Pangenome with gene families
    :param source: source of the metadata
    :param metatype: select to which pangenome element metadata will be added
    :param omit: allow to omit a row in dataframe if the element name is not find in pangenomes
    :param disable_bar: Disable progress bar

    :raise KeyError: element name is not find in pangenome
    :raise AssertionError: Metatype is not recognized
    """
    assert metatype in ["families", "genomes", "genes", "RGPs", "spots", "modules"]
    for row in tqdm(metadata_df.iterrows(), unit='row',
                    total=metadata_df.shape[0], disable=disable_bar):
        row = row[1]
        try:
            if metatype == "families":
                element = pangenome.get_gene_family(row[metatype])
            elif metatype == "genomes":
                element = pangenome.get_organism(row[metatype])
            elif metatype == "genes":
                element = pangenome.get_gene(row[metatype])
            elif metatype == "RGPs":
                element = pangenome.get_region(row[metatype])
            elif metatype == "spots":
                element = pangenome.get_spot(row[metatype])
            else:  # metatype == "modules":
                element = pangenome.get_module(row[metatype])
        except KeyError:
            if not omit:
                raise KeyError(f"{metatype} {row[metatype]} does not exist in pangenome. Check name in your file")
            else:
                logging.getLogger().debug(f"{metatype}: {row[metatype]} doesn't exist")
        else:
            meta = Metadata(source=source, **{k: v for k, v in row.to_dict().items() if k != metatype})
            element.add_metadata(source=source, metadata=meta)

    pangenome.status["metadata"][metatype] = "Computed"
    pangenome.status["metasources"][metatype].append(source)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    metadata_df = check_metadata_format(args.metadata, args.assign)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    check_pangenome_metadata(pangenome, source=args.source, metatype=args.assign,
                             force=args.force, disable_bar=args.disable_prog_bar)
    assign_metadata(metadata_df, pangenome=pangenome, source=args.source, metatype=args.assign,
                    omit=args.omit, disable_bar=args.disable_prog_bar)
    logging.getLogger().info("Metadata assignment Done")
    write_pangenome(pangenome, pangenome.file, force=args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("metadata", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_meta(parser)
    return parser


def parser_meta(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of graph command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-m', '--metadata', type=Path, nargs='?',
                          help='Metadata in TSV file. See our github for more detail about format')
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the metadata source')
    required.add_argument("-a", "--assign", required=True, type=str, nargs="?",
                          choices=["families", "genomes", "genes", "RGPs", "spots", "modules"],
                          help="Select to which pangenome element metadata will be assigned")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--omit", required=False, action="store_true",
                          help="Allow to pass if a key in metadata is not find in pangenome")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_meta(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
