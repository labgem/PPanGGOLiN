#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
from pathlib import Path
import csv
from typing import Union

# installed libraries
from tqdm import tqdm
import pandas as pd

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.metadata import Metadata
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome


def check_pangenome_metadata(pangenome: Pangenome, source: str, metatype: str, force: bool = False,
                             disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param source: source of the annotation
    :param force: erase if an annotation for the provide source already exist
    :param disable_bar: Disable bar
    """
    need_annotations = False
    need_families = False
    if metatype in ["genomes", "genes", "families"]:
        need_annotations = True
    if metatype == "families":
        need_families = True

    if pangenome.status["metadata"][metatype] == "inFile" and source in pangenome.status["metasources"][metatype]:
        if force:
            erase_pangenome(pangenome, metadata=True, source=source, metatype=metatype)
        else:
            raise Exception(
                f"An annotation corresponding to the source : '{source}' already exist in pangenome organims."
                "Add the option --force to erase")
    check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                         disable_bar=disable_bar)


def assign_metadata_to_families(metadata: Path, pangenome: Pangenome, source: str,
                                omit: bool = False, disable_bar: bool = False):
    meta_col_names = ['gene_family', 'protein_name', 'accession', 'e_value',
                      'score', 'bias', 'secondary_name', 'description']
    meta_col_type = {'gene_family': str,
                     'protein_name': str,
                     'accession': str,
                     'e_value': float,
                     'score': float,
                     'bias': float,
                     'secondary_name': str,
                     'description': str}
    metadata_df = pd.read_csv(metadata, sep="\t", header=None, quoting=csv.QUOTE_NONE,
                              names=meta_col_names, dtype=meta_col_type)
    for row in tqdm(metadata_df.itertuples(index=False), unit='row',
                    total=metadata_df.shape[0], disable=disable_bar):
        try:
            gene_families = pangenome.get_gene_family(name=row.gene_families)
        except KeyError:
            if not omit:
                raise KeyError(f"Family {row.Gene_family} does not exist in pangenome. Check name in your file")
        else:
            annotation = Metadata(source=source, accession=row.accession, value=row.value,
                                  description=row.description, score=row.score, e_val=row.e_value, bias=row.bias)
            gene_families.add_metadata(source=source, metadata=annotation)

    pangenome.status["metadata"]["families"] = "Computed"
    pangenome.status["metasources"]["families"].append(source)


def assign_metadata_to_genomes(metadata: Path, pangenome: Pangenome, source: str,
                               omit: bool = False, disable_bar: bool = False):
    meta_col_type = {'genome': str,
                     'value': str,
                     'accession': str,
                     'e_value': float,
                     'score': float,
                     'bias': float,
                     'Description': str}
    metadata_df = pd.read_csv(metadata, sep="\t", header=0, quoting=csv.QUOTE_NONE, dtype=meta_col_type)
    if not {"genome", "value"}.issubset(set(metadata_df.columns)):
        raise KeyError("You should at least provide in columns names : genome and value."
                       "Look at documentation for more information")
    elif not set(metadata_df.columns).issubset(set(meta_col_type.keys())):
        raise KeyError("There are one or more column name not acceptable."
                       f"Acceptable names are : {list(meta_col_type.keys())}")
    pd.to_numeric(metadata_df["value"], downcast='integer', errors='ignore')
    for row in tqdm(metadata_df.itertuples(index=False), unit='row',
                    total=metadata_df.shape[0], disable=disable_bar):
        try:
            genome = pangenome.get_organism(row.genome)
        except KeyError:
            if not omit:
                raise KeyError(f"Genome {row.genome} does not exist in pangenome. Check name in your file")
        else:
            meta = Metadata(source=source, value=row.value,
                                  accession=row.accession if "accession" in metadata_df.columns else None,
                                  description=row.description if "description" in metadata_df.columns else None,
                                  score=row.score if "score" in metadata_df.columns else None,
                                  e_val=row.score if "e_val" in metadata_df.columns else None,
                                  bias=row.bias if "bias" in metadata_df.columns else None)
            genome.add_metadata(source=source, metadata=meta)

    pangenome.status["metadata"]["genomes"] = "Computed"
    pangenome.status["metasources"]["genomes"].append(source)


def assign_metadata_to_genes(metadata: Path, pangenome: Pangenome, source: str,
                             omit: bool = False, disable_bar: bool = False):
    meta_col_names = ['genes', 'value', 'accession', 'e_value',
                      'score', 'bias', 'description']
    meta_col_type = {'genes': str,
                     'value': Union[str, int, float],
                     'accession': str,
                     'e_value': float,
                     'score': float,
                     'bias': float,
                     'Description': str}
    metadata_df = pd.read_csv(metadata, sep="\t", header=None, quoting=csv.QUOTE_NONE,
                              names=meta_col_names, dtype=meta_col_type)
    for row in tqdm(metadata_df.itertuples(index=False), unit='row',
                    total=metadata_df.shape[0], disable=disable_bar):
        try:
            gene = pangenome.get_gene(row.genes)
        except KeyError:
            if not omit:
                raise KeyError(f"Gene {row.gene} does not exist in pangenome. Check name in your file")
        else:
            annotation = Metadata(source=source, accession=row.accession, value=row.value,
                                  description=row.description, score=row.score, e_val=row.e_value, bias=row.bias)
            gene.add_metadata(source=source, metadata=annotation)

    pangenome.status["metadata"]["genomes"] = "Computed"
    pangenome.status["metasources"]["genomes"].append(source)


def assign_metadata(metadata: Path, pangenome: Pangenome, source: str, metatype: str,
                    omit: bool = False, disable_bar: bool = False) -> dict:
    """ Add to gene families an annotation and create a dictionary with for each annotation a set of gene family

    :param metadata: Dataframe with for each family an annotation
    :param pangenome: Pangenome with gene families
    :param source: source of the annotation
    :param disable_bar:
    :return: Dictionary with for each annotation a set of gene family
    """
    if metatype == "families":
        assign_metadata_to_families(metadata, pangenome, source, omit, disable_bar)
    elif metatype == "genomes":
        assign_metadata_to_genomes(metadata, pangenome, source, omit, disable_bar)
    elif metatype == "genes":
        assign_metadata_to_genes(metadata, pangenome, source, omit, disable_bar)
    elif metatype == "RGPs":
        raise NotImplementedError("Option not implemented yet !")
    elif metatype == "spots":
        raise NotImplementedError("Option not implemented yet !")
    elif metatype == "modules":
        raise NotImplementedError("Option not implemented yet !")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    check_pangenome_metadata(pangenome, source=args.source, metatype=args.assign,
                             force=args.force, disable_bar=args.disable_prog_bar)
    assign_metadata(metadata=args.metadata, pangenome=pangenome, source=args.source, metatype=args.assign,
                    omit=args.omit, disable_bar=args.disable_prog_bar)
    logging.getLogger().info("Metadata assignment Done")
    write_pangenome(pangenome, pangenome.file, disable_bar=args.disable_prog_bar)


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
                          help='Gene families annotation in TSV file. See our github for more detail about format')
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source. Default use name of annnotation file or directory.')
    required.add_argument("-a", "--assign", required=True, type=str, nargs="?",
                          choices=["families", "genomes", "genes", "RGPs", "spots", "modules"],
                          help="Select to which pangenome element metadata will be assigned.")
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
