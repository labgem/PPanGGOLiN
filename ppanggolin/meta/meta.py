#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
from pathlib import Path
import csv

# installed libraries
from tqdm import tqdm
import pandas as pd

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.metadata import Metadata
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome

res_col_names = ['Gene_family', 'Accession', 'protein_name', 'e_value',
                 'score', 'bias', 'secondary_name', 'Description']


def check_pangenome_annotation(pangenome: Pangenome, source: str, force: bool = False, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param source: source of the annotation
    :param force: erase if an annotation for the provide source already exist
    :param disable_bar: Disable bar
    """
    if "annotations_sources" in pangenome.status and source in pangenome.status["annotations_sources"]:
        if force:
            erase_pangenome(pangenome, annotations=True, source=source)
        else:
            raise Exception(f"An annotation corresponding to the source : '{source}' already exist in pangenome."
                            f"Add the option --force to erase")
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=disable_bar)


def annotation_to_families(annotation_df: pd.DataFrame, pangenome: Pangenome, source: str = None) -> dict:
    """ Add to gene families an annotation and create a dictionary with for each annotation a set of gene family

    :param annotation_df: Dataframe with for each family an annotation
    :param pangenome: Pangenome with gene families
    :param source: source of the annotation

    :return: Dictionary with for each annotation a set of gene family
    """
    for row in annotation_df.itertuples(index=False):
        gene_fam = pangenome.get_gene_family(name=row.Gene_family)
        if gene_fam is not None:
            annotation = Metadata(source=source, accession=row.Accession, value=row.protein_name,
                                  secondary_names=row.secondary_name, description=row.Description,
                                  score=row.score, e_val=row.e_value, bias=row.bias)
            gene_fam.add_metadata(source=source, metadata=annotation)
        else:
            logging.getLogger().warning(f"Family {row.Gene_family} does not exist in pangenome. "
                                        "Check name in your file")
        pangenome.status["annotations"] = "Computed"


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    check_pangenome_annotation(pangenome, source=args.source, force=args.force, disable_bar=args.disable_prog_bar)
    annotation_df = pd.read_csv(args.data, sep="\t", header=None, quoting=csv.QUOTE_NONE, names=res_col_names)
    annotation_to_families(annotation_df, pangenome, args.source)
    logging.getLogger().info("Annotation Done")
    logging.getLogger().info("Write Annotation in pangenome")
    # write_pangenome(pangenome, pangenome_info["path"], source=args.source, disable_bar=args.disable_prog_bar)


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
    required.add_argument('-d', '--data', type=Path, nargs='?',
                          help='Gene families annotation in TSV file. See our github for more detail about format')
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source. Default use name of annnotation file or directory.')
    # optional = parser.add_argument_group(title="Optional arguments")


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
