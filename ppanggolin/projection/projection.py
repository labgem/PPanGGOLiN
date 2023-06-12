#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
from multiprocessing import get_context
import logging
import os
import time

# installed libraries
from tqdm import tqdm

# # local libraries
# from ppanggolin.annotate.synta import annotate_organism, read_fasta, get_dna_sequence
# from ppanggolin.pangenome import Pangenome
# from ppanggolin.genome import Organism, Gene, RNA, Contig
# from ppanggolin.utils import read_compressed_or_not, mk_file_name, min_one, restricted_float
# from ppanggolin.formats import write_pangenome


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pass


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("projection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_projection(parser)
    return parser


def parser_projection(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of annotate command

    :param parser: parser for annotate argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=False, type=str, help="The pangenome.h5 file")
    
    required.add_argument('--organism_name', required=False, type=str,
                        help="Name of the organism whose genome is being projected onto the provided pangenome.")
    
    required.add_argument('--fasta_file', required=False, type=str,
                        help="The filepath of the genomic sequence(s) in FASTA format for the projected genome. "
                        "(Fasta file can be compressed with gzip)")

    required.add_argument('--anno_file', required=False, type=str,
                        help="The filepath of the annotations in GFF/GBFF format for the projected genome. "
                        "(Annotation file can be compressed with gzip)")

    # required.add_argument('--fasta', required=False, type=str,
    #                       help="A tab-separated file listing the organism names, and the fasta filepath of its genomic "
    #                            "sequence(s) (the fastas can be compressed with gzip). One line per organism.")
    
    # required.add_argument('--anno', required=False, type=str,
    #                       help="A tab-separated file listing the organism names, and the gff/gbff filepath of its "
    #                            "annotations (the files can be compressed with gzip). One line per organism. "
    #                            "If this is provided, those annotations will be used.")

    annotate = parser.add_argument_group(title="Annotation arguments")

    annotate.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    annotate.add_argument('--allow_overlap', required=False, action='store_true', default=False,
                          help="Use to not remove genes overlapping with RNA features.")
    annotate.add_argument("--norna", required=False, action="store_true", default=False,
                          help="Use to avoid annotating RNA features.")
    annotate.add_argument("--kingdom", required=False, type=str.lower, default="bacteria",
                          choices=["bacteria", "archaea"],
                          help="Kingdom to which the prokaryota belongs to, "
                               "to know which models to use for rRNA annotation.")
    annotate.add_argument("--translation_table", required=False, type=int, default=11,
                          help="Translation table (genetic code) to use.")
    annotate.add_argument("--basename", required=False, default="pangenome", help="basename for the output file")

    annotate.add_argument("--prodigal_procedure", required=False, type=str.lower, choices=["single", "meta"],
                          default=None, help="Allow to force the prodigal procedure. "
                                             "If nothing given, PPanGGOLiN will decide in function of contig length")
    

    cluster = parser.add_argument_group(title="Clustering arguments")
    cluster.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family.")
    cluster.add_argument('--identity', required=False, type=float, default=0.5,
                          help="min identity percentage threshold")
    cluster.add_argument('--coverage', required=False, type=float, default=0.8,
                          help="min coverage percentage threshold")
    cluster.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    
