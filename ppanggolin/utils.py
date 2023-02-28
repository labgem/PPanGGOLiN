#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys
import os
import gzip
import argparse
from io import TextIOWrapper
from pathlib import Path
from typing import TextIO, Union, BinaryIO

import networkx as nx
import pkg_resources
from numpy import repeat
import logging
import tempfile

from scipy.sparse import csc_matrix

import yaml
from collections import defaultdict


from ppanggolin.geneFamily import GeneFamily


def check_log(name: str) -> TextIO:
    """Check if the output log is writable

    :param name: Path to the log output

    :return: output for log
    """
    if name == "stdout":
        return sys.stdout
    elif name == "stderr":
        return sys.stderr
    else:
        try:
            log_file = open(name, "w")
        except IOError:
            raise IOError("The given log file does not appear.")
        except Exception:
            raise Exception("An unexpected error happened with your logfile. Please check if he is accessible."
                            "If everything looks good, please report an issue on our GitHub.")
        else:
            return log_file


def check_tsv_sanity(tsv):
    """ Check if the given tsv is readable for the next PPanGGOLiN step

    :param tsv: Path to the input tsv
    """
    f = open(tsv, "r")
    name_set = set()
    duplicated_names = set()
    non_existing_files = set()
    for line in f:
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:
            raise Exception(f"No tabulation separator found in given file: {tsv}")
        if " " in elements[0]:
            raise Exception(f"Your genome names contain spaces (The first encountered genome name that had this string:"
                            f" '{elements[0]}'). To ensure compatibility with all of the dependencies of PPanGGOLiN "
                            f"this is not allowed. Please remove spaces from your genome names.")
        old_len = len(name_set)
        name_set.add(elements[0])
        if len(name_set) == old_len:
            duplicated_names.add(elements[0])
        if not os.path.exists(elements[1]):
            non_existing_files.add(elements[1])
    if len(non_existing_files) != 0:
        raise Exception(f"Some of the given files do not exist. The non-existing files are the following : "
                        f"'{' '.join(non_existing_files)}'")
    if len(duplicated_names) != 0:
        raise Exception(f"Some of your genomes have identical names. The duplicated names are the following : "
                        f"'{' '.join(duplicated_names)}'")


def check_input_files(anno: str = None, pangenome: str = None, fasta: str = None):
    """ Checks if the provided input files exist and are of the proper format

    :param anno: Path to the annotation file
    :param pangenome: Path to the pangenome hdf5 file
    :param fasta: path to the fasta file
    """
    if pangenome is not None and not os.path.exists(pangenome):
        raise FileNotFoundError(f"No such file or directory: '{pangenome}'")

    if anno is not None:
        if not os.path.exists(anno):
            raise FileNotFoundError(f"No such file or directory: '{anno}'")
        check_tsv_sanity(anno)

    if fasta is not None:
        if not os.path.exists(fasta):
            raise FileNotFoundError(f"No such file or directory: '{fasta}'")
        check_tsv_sanity(fasta)


def set_verbosity_level(args):
    """Set the verbosity level

    :param args: argument pass by command line
    """
    level = logging.INFO  # info, warnings and errors, default verbose == 1
    if hasattr(args, "verbose"):
        if args.verbose == 2:
            level = logging.DEBUG  # info, debug, warnings and errors
        elif args.verbose == 0:
            level = logging.WARNING  # only warnings and errors

        if args.log != sys.stdout and not args.disable_prog_bar:  # if output is not to stdout we remove progress bars.
            args.disable_prog_bar = True

        logging.basicConfig(stream=args.log, level=level,
                            format='%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
        logging.getLogger().info("Command: " + " ".join([arg for arg in sys.argv]))
        logging.getLogger().info("PPanGGOLiN version: " + pkg_resources.get_distribution("ppanggolin").version)


def jaccard_similarities(mat: csc_matrix, jaccard_similarity_th) -> csc_matrix:
    """ Compute the jaccard similarities

    :param mat:
    :param jaccard_similarity_th: threshold

    :return:
    """
    cols_sum = mat.getnnz(axis=0)
    ab = mat.T * mat
    # for rows
    aa = repeat(cols_sum, ab.getnnz(axis=0))
    # for columns
    bb = cols_sum[ab.indices]
    similarities = ab.copy()
    similarities.data /= (aa + bb - ab.data)
    similarities.data[similarities.data < jaccard_similarity_th] = 0
    similarities.eliminate_zeros()
    return similarities


def read_compressed_or_not(file_or_file_path: Union[str, BinaryIO, TextIOWrapper, TextIO]) -> Union[TextIOWrapper,
                                                                                                    BinaryIO, TextIO]:
    """
    Reads a file object or file path, uncompresses it, if need be.

    :param file_or_file_path: Path to the input file

    :return: TextIO object in read only
    """
    file = file_or_file_path
    if isinstance(file, str):
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except AttributeError:
            return file
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        return TextIOWrapper(gzip.open(filename=file, mode="r"))
    else:
        file.close()
        file = open(file.name, "r")
        return file


def write_compressed_or_not(file_path: str, compress: bool = False) -> Union[gzip.GzipFile, TextIO]:
    """
    Create a file-like object, compressed or not.

    :param file_path: Path to the file
    :param compress: Compress the file in .gz

    :return: file-like object, compressed or not
    """
    if compress:
        return gzip.open(file_path + ".gz", mode="wt")
    else:
        return open(file_path, "w")


def is_compressed(file_or_file_path: Union[str, TextIO, gzip.GzipFile]):
    """ Checks is a file, or file path given is compressed or not

    :param file_or_file_path: Input file

    :return: Get if the file is compressed
    """
    file = file_or_file_path
    if isinstance(file, str):
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except AttributeError:
            return False
    if file.read(2).startswith(b'\x1f\x8b'):
        return True
    file.close()
    return False


def mk_outdir(output, force):
    """ Create a directory at the given output if it doesn't exist already

    :param output: Path where to create directory
    :param force: Force to write in the directory

    :raise FileExistError: The current path already exist and force is false
    """
    if not os.path.exists(output):
        os.makedirs(output)
    elif not force:
        raise FileExistsError(f"{output} already exists. Use -f if you want to overwrite the files in the directory")


def mk_file_name(basename: str, output: str, force: bool = False) -> Path:
    """Returns a usable filename for a ppanggolin output file, or crashes.

    :param basename: basename for the file
    :param output: Path to save the file
    :param force: Force to write the file

    :return: Path to the file
    """
    filename = Path(output + "/" + basename)
    if filename.suffix != ".h5":
        filename = filename.with_suffix(".h5")

    mk_outdir(output, force)

    if filename.exists() and not force:
        raise FileExistsError(f"{filename.name} already exists. Use -f if you want to overwrite the file")
    return filename


def restricted_float(x) -> float:
    """Decrease the choice possibility of float in argparse

    :param x: given float by user

    :return: given float if it is acceptable

    :raise argparse.ArgumentTypeError: The float is not acceptable
    """
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def min_one(x) -> int:
    """Check if the given int is superior to one

    :param x: given float by user

    :return: given float if it is acceptable

    :raise argparse.ArgumentTypeError: The float is not acceptable
    """
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("%r is inferior to 1" % (x,))
    return x


def connected_components(g: nx.Graph, removed: set, weight: float):
    """
    Yields subgraphs of each connected component you get when filtering edges based on the given weight.

    :param g: Subgraph
    :param removed: removed node
    :param weight: threshold to remove node or not
    """
    for v in g.nodes:
        if v not in removed:
            c = set(_plain_bfs(g, v, removed, weight))
            yield c
            removed.update(c)


def _plain_bfs(g: nx.Graph, source: GeneFamily, removed: set, weight: float):
    """
    A fast BFS node generator, copied from networkx then adapted to the current use case

    :param g: graph with the nodes
    :param source: current node
    :param removed: set of removed nodes
    :param weight:threshold to remove node or not
    """

    nextlevel = {source}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in removed:
                yield v
                removed.add(v)

                for n in g.neighbors(v):
                    if n not in removed:
                        edge_genes_v = g[v][n]["genes"][v]
                        edge_genes_n = g[v][n]["genes"][n]
                        # if the edge is indeed existent for most genes of both families, we use it
                        if len(edge_genes_n) / len(g.nodes[n]["genes"]) >= weight and len(edge_genes_v) / len(
                                g.nodes[v]["genes"]) >= weight:
                            nextlevel.add(n)


def add_gene(obj, gene, fam_split: bool = True):
    """

    :param obj:
    :param gene:
    :param fam_split:
    """
    if fam_split:
        try:
            obj["genes"][gene.family].add(gene)
        except KeyError:
            try:
                obj["genes"][gene.family] = {gene}
            except KeyError:
                obj["genes"] = {gene.family: {gene}}
    else:
        try:
            obj["genes"].add(gene)
        except KeyError:
            obj["genes"] = {gene}


def check_option_workflow(args):
    """
    Check if the given argument to a workflow command is usable

    :param args: list of arguments
    """
    if args.clusters is not None and not any([args.fasta, args.anno]):
        raise Exception("If you give --clusters option, you must give at least --fasta or --anno")

    if not any([args.fasta, args.anno]):
        raise Exception("At least one of --fasta or --anno must be given")


def parse_config_file(yaml_config_file: str) -> dict:
    """
    Parse yaml config file.

    :param yaml_config_file: config file in yaml

    :return: dict of config with key the command and as value another dict with param as key and value as value. 
    """

    with yaml_config_file as yaml_fh:
        config = yaml.safe_load(yaml_fh)
        
    return config

def overwrite_params_with_cli_args(config_args: argparse.Namespace, cli_args: argparse.Namespace):
    """
    Overwrite arguments from config/default with cli arguments.
    When arguments are given in CLI, their value is used instead of the one found in config or the default one.
    CLI argument that have not been given in CLI have a None value. 

    :param config_args: Arguments with config values or default values.
    :param cli_args: Arguments parsed from the cmd line.

    :return: object with arguments
    """

    non_default_args = [arg for arg in dir(cli_args) if not arg.startswith('_')]

    for non_default_arg in non_default_args:
        arg_val = getattr(cli_args, non_default_arg)
        
        if non_default_arg in config_args and arg_val:
            logging.getLogger().info(f'Parameter "--{non_default_arg} {arg_val}" has been specified in command line. Its value overwrite putative config value.')
            setattr(config_args, non_default_arg, arg_val)
    

def get_cmd_args_from_config(step_name: str, parser_fct, config_param_val: dict,
                            cli_args: argparse.Namespace,
                            general_params: list) -> argparse.Namespace: 
    """
    Parse arguments from config file of a specific command using argument parser of this command. 

    :param step_name: name of the step (ie annotate, cluster.. )
    :param parser_fct: parser function of the command
    :param config_param_val: dict parsed from config file with key value for the command
    :param cli_args: Arguments parsed from the cmd line which overwrite config or default values when they are specified in cmd line.
    :param general_params: General parameters to remove from the expected arguments. These parameters are managed by cmd line arguments directly.

    :return: object with arguments for the given command 
    """

    # Abbreviations are not allowed in config file
    parser = argparse.ArgumentParser(prog=f"{step_name} args from config file", 
                                    allow_abbrev=False, add_help=False)  

    parser_fct(parser)

    # remove required arguments. Config file can only contained optional arguments
    parser._actions = [p_action for p_action in parser._actions if p_action.required == False]

    # remove general arguments to only expect arguments specific to the step.
    parser._actions = [p_action for p_action in parser._actions if p_action.dest not in general_params]
    
    parser.usage = "Yaml expected structure"
    
    arguments_to_parse = []
    off_flags = []
    for param, val in config_param_val.items():

        if type(val) == bool:
            # param is a flag
            if val is True:
                arguments_to_parse.append(f"--{param}")
            else:
                off_flags.append(param)
            # if val is false, the param is not added to the list
        else:
            # argument is a "--param val" type
            arguments_to_parse.append(f"--{param}")

            if type(val) == list:
                # range of values need to be added one by one
                arguments_to_parse += [str(v) for v in val]
            else:
                arguments_to_parse.append(str(val))

    args = parser.parse_args(arguments_to_parse)

    logging.getLogger().info(f'Config file: {step_name}: {len(config_param_val)} arguments parsed from config.')

    if config_param_val:
        logging.getLogger().debug(f'Arguments to parse: {" ".join(arguments_to_parse)}')
        logging.getLogger().debug(f'Flag arguments set to False: {off_flags}')

    overwrite_params_with_cli_args(args, cli_args)

    return args


def add_common_arguments(subparser: argparse.ArgumentParser):
    """
    Add common argument to the input subparser.

    :param subparser: A subparser object from any subcommand.
    """

    common = subparser._action_groups.pop(1)  # get the 'optional arguments' action group.
    common.title = "Common arguments"
    common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(),
                        help="directory for storing temporary files")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    # common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")

    subparser._action_groups.append(common)
