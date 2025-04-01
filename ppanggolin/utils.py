#!/usr/bin/env python3

# default libraries
import logging
import sys
import os
import gzip
import bz2
import zipfile
import argparse
import inspect

from io import TextIOWrapper
from pathlib import Path
from typing import (
    Collection,
    TextIO,
    Union,
    BinaryIO,
    Tuple,
    List,
    Set,
    Iterable,
    Dict,
    Any,
    Optional,
)
from contextlib import contextmanager
import tempfile
import time
from itertools import zip_longest
import re
import subprocess
import shutil

import networkx as nx
from importlib.metadata import distribution
from numpy import repeat
from collections.abc import Callable

from scipy.sparse import csc_matrix

import yaml
from collections import defaultdict

# all input params that exists in ppanggolin
ALL_INPUT_PARAMS = [
    "fasta",
    "anno",
    "clusters",
    "pangenome",
    "fasta_file",
    "annot_file",
    "genome_name",
]  # the last three params is for projection cmd

# all params that should be in the general_parameters section of the config file
ALL_GENERAL_PARAMS = [
    "output",
    "basename",
    "rarefaction",
    "no_flat_files",
    "tmpdir",
    "verbose",
    "log",
    "disable_prog_bar",
    "force",
    "config",
]

WORKFLOW_SUBCOMMANDS = {"all", "workflow", "panrgp", "panmodule"}

# command that can be launched inside a workflow subcommand
ALL_WORKFLOW_DEPENDENCIES = [
    "annotate",
    "cluster",
    "graph",
    "partition",
    "rarefaction",
    "rgp",
    "spot",
    "module",
    "draw",
    "write_pangenome",
    "write_genomes",
]

# Inside a workflow command, write output default is overwrite to output some flat files
WRITE_PAN_FLAG_DEFAULT_IN_WF = [
    "csv",
    "Rtab",
    "gexf",
    "light_gexf",
    "stats",
    "json",
    "partitions",
    "regions",
    "borders",
    "modules",
    "spot_modules",
    "spots",
    "families_tsv",
]
WRITE_GENOME_FLAG_DEFAULT_IN_WF = ["table", "proksee", "gff"]

DRAW_FLAG_DEFAULT_IN_WF = ["tile_plot", "ucurve", "draw_spots"]


def check_log(log_file: str) -> TextIO:
    """
    Check if the output log is writable

    :param log_file: Path to the log output

    :return: output for log
    """
    if log_file == "stdout":
        return sys.stdout
    elif log_file == "stderr":
        return sys.stderr

    elif os.path.exists(log_file):
        # path exists
        if os.path.isfile(log_file):  # is it a file or a dir?
            # also works when file is a link and the target is writable
            if os.access(log_file, os.W_OK):
                return log_file
            else:
                raise OSError(
                    f"The given log file {log_file} is not writable. Please check if it is accessible."
                )
        else:
            raise OSError(
                f"The given log file: {log_file} is a directory. Please provide a valid log file."
            )

    # target does not exist, check perms on parent dir
    parent_dir = os.path.dirname(log_file)
    if not parent_dir:
        parent_dir = "."
    # target is creatable if parent dir is writable
    if os.access(parent_dir, os.W_OK):
        return log_file
    else:
        raise OSError(
            f"The given log file {log_file} is not writable. Please check if it is accessible."
        )


def check_tsv_sanity(tsv_file: Path):
    """Check if the given TSV file is readable for the next PPanGGOLiN step.

    :param tsv: Path to the TSV containing organism information.
    :raises ValueError: If the file format is incorrect or contains invalid genome names.
    """
    with read_compressed_or_not(tsv_file) as input_file:
        name_set = set()
        duplicated_names = set()
        non_existing_files = set()

        for line in input_file:
            elements = [el.strip() for el in line.split("\t")]

            if len(elements) <= 1:
                raise ValueError(f"No tabulation separator found in file: {tsv_file}")

            genome_name, genome_path = elements[0], elements[1]

            if " " in genome_name:
                raise ValueError(
                    f"Genome names cannot contain spaces (first encountered: '{genome_name}'). "
                    "Please remove spaces to ensure compatibility with PPanGGOLiN dependencies."
                )

            if genome_name in name_set:
                duplicated_names.add(genome_name)
            name_set.add(genome_name)

            org_path = Path(genome_path)
            if (
                not org_path.exists()
                and not tsv_file.parent.joinpath(org_path).exists()
            ):
                non_existing_files.add(genome_path)

        if non_existing_files:
            raise ValueError(
                f"Some specified genome files do not exist: {', '.join(non_existing_files)}"
            )
        if duplicated_names:
            raise ValueError(
                f"Some genome names are duplicated: {', '.join(duplicated_names)}"
            )


def check_input_files(file: Path, check_tsv: bool = False):
    """Checks if the provided input files exist and are of the proper format

    :param file: Path to the file
    :param check_tsv: Allow checking tsv file for annotation or fasta list
    """
    if file.exists():
        if check_tsv:
            check_tsv_sanity(file)
    else:
        raise FileNotFoundError(
            f"No such file or directory: '{file.absolute().as_posix()}'"
        )


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

        if (
            args.log != sys.stdout and not args.disable_prog_bar
        ):  # if output is not to stdout we remove progress bars.
            args.disable_prog_bar = True
        str_format = "%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s"
        datefmt = "%Y-%m-%d %H:%M:%S"
        if args.log in [sys.stdout, sys.stderr]:
            # use stream
            logging.basicConfig(
                stream=args.log, level=level, format=str_format, datefmt=datefmt
            )
        else:
            # log is written in a files. basic condif uses filename
            logging.basicConfig(
                filename=args.log, level=level, format=str_format, datefmt=datefmt
            )
        logging.getLogger("PPanGGOLiN").info(
            "Command: " + " ".join(arg for arg in sys.argv)
        )
        logging.getLogger("PPanGGOLiN").info(
            f"PPanGGOLiN version: {distribution('ppanggolin').version}"
        )


def jaccard_similarities(mat: csc_matrix, jaccard_similarity_th) -> csc_matrix:
    """Compute the jaccard similarities

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
    similarities.data /= aa + bb - ab.data
    similarities.data[similarities.data < jaccard_similarity_th] = 0
    similarities.eliminate_zeros()
    return similarities


def is_compressed(
    file_or_file_path: Union[Path, BinaryIO, TextIOWrapper, TextIO],
) -> Tuple[bool, Union[str, None]]:
    """
    Detects if a file is compressed based on its file signature.

    :param file_or_file_path: The file to check.

    :return: True if the file is a recognized compressed format with the format name, False otherwise.

    :raises TypeError: If the file type is not supported.
    """
    file_signatures = {
        b"\x1f\x8b": "gzip",
        b"BZh": "bz2",
        b"\x50\x4b\x03\x04": "zip",
        b"\xfd\x37\x7a\x58\x5a\x00": "xz",
    }

    def check_file_signature(byte_stream) -> Tuple[bool, Union[str, None]]:
        """
        Checks if the provided byte stream starts with a known file signature.

        :param byte_stream: The first few bytes of a file.

        :return: True if the byte stream starts with a known file signature, False otherwise.
        """
        for signature, filetype in file_signatures.items():
            if byte_stream.startswith(signature):
                return True, filetype
        return False, None

    # Determine the type of file and read its first few bytes
    if isinstance(file_or_file_path, Path):
        with file_or_file_path.open("rb") as file:
            first_bytes = file.read(4)
    else:
        if isinstance(file_or_file_path, BinaryIO):
            first_bytes = file_or_file_path.readline()[:4]
        elif isinstance(file_or_file_path, TextIOWrapper):
            first_bytes = file_or_file_path.buffer.read(4)
        elif isinstance(file_or_file_path, TextIO):
            first_bytes = file_or_file_path.read(4).encode()
        else:
            raise TypeError("Unsupported file type")
        file_or_file_path.seek(0)  # Reset the file position

    return check_file_signature(first_bytes)


def read_compressed_or_not(
    file_or_file_path: Union[Path, BinaryIO, TextIOWrapper, TextIO],
) -> Union[TextIOWrapper, BinaryIO, TextIO]:
    """
    Opens and reads a file, decompressing it if necessary.

    Parameters:
    file (pathlib.Path, io.BytesIO, io.TextIOWrapper, io.TextIOBase): The file to read.
    It can be a Path object from the pathlib module, a BytesIO object, a TextIOWrapper, or TextIOBase object.

    Returns:
    str: The contents of the file, decompressed if it was a recognized compressed format.

    Raises:
    TypeError: If the file type is not supported.
    """
    is_comp, comp_type = is_compressed(file_or_file_path)
    if is_comp:
        if comp_type == "gzip":
            return gzip.open(file_or_file_path, "rt")
        elif comp_type == "bz2":
            return bz2.open(file_or_file_path, "rt")
        elif comp_type == "xz":
            raise NotImplementedError(
                "Unfortunately PPanGGOLiN does not support xz compressed files. "
                "Please report an issue on our GitHub to let us know we should work on it."
            )
        elif comp_type == "zip":
            with zipfile.ZipFile(file_or_file_path, "r") as z:
                logging.getLogger("PPanGGOLiN").warning(
                    "Assuming we want to read the first file in the ZIP archive"
                )
                file_list = z.namelist()
                if file_list:
                    return TextIOWrapper(z.open(file_list[0], "r"))
    else:  # Non-compressed file
        if isinstance(file_or_file_path, Path):
            return open(file_or_file_path)
        else:
            return file_or_file_path


def write_compressed_or_not(
    file_path: Path, compress: bool = False
) -> Union[gzip.GzipFile, TextIOWrapper]:
    """
    Create a file-like object, compressed or not.

    :param file_path: Path to the file
    :param compress: Compress the file in .gz

    :return: file-like object, compressed or not
    """
    if compress:
        return gzip.open(file_path.parent / (file_path.name + ".gz"), mode="wt")
    else:
        return open(file_path, "w")


def mk_outdir(output: Path, force: bool = False, exist_ok: bool = False):
    """Create a directory at the given output if it doesn't exist already

    :param output: Path where to create directory
    :param force: Force to write in the directory
    :param exist_ok: Does not give an error if the directory already exists.

    :raise FileExistError: The current path already exist and force is false
    """
    if not output.is_dir():
        logging.getLogger("PPanGGOLiN").debug(
            f"Create output directory {output.absolute().as_posix()}"
        )
        Path.mkdir(output, exist_ok=exist_ok)
    else:
        if not force:
            raise FileExistsError(
                f"{output} already exists. Use -f if you want to overwrite the files in the directory"
            )


@contextmanager
def create_tmpdir(main_dir, basename="tmpdir", keep_tmp=False):
    if keep_tmp:
        dir_name = basename + time.strftime("_%Y-%m-%d_%H.%M.%S", time.localtime())

        new_tmpdir = main_dir / dir_name
        logging.getLogger("PPanGGOLiN").debug(
            f"Creating a temporary directory: {new_tmpdir.as_posix()}. This directory will be retained."
        )

        mk_outdir(new_tmpdir, force=True)
        yield new_tmpdir

    else:
        with tempfile.TemporaryDirectory(dir=main_dir, prefix=basename) as new_tmpdir:
            logging.getLogger("PPanGGOLiN").debug(
                f"Creating a temporary directory: {new_tmpdir}. This directory won't be retained."
            )
            yield Path(new_tmpdir)


def mk_file_name(basename: str, output: Path, force: bool = False) -> Path:
    """Returns a usable filename for a ppanggolin output file, or crashes.

    :param basename: basename for the file
    :param output: Path to save the file
    :param force: Force to write the file

    :return: Path to the file
    """
    filename = output / basename
    if filename.suffix != ".h5":
        filename = filename.with_suffix(".h5")

    mk_outdir(output, force)

    if filename.exists() and not force:
        raise FileExistsError(
            f"{filename.name} already exists. Use -f if you want to overwrite the file"
        )
    return filename


def detect_filetype(filename: Path) -> str:
    """
    Detects whether the current file is gff3, gbk/gbff, fasta, tsv or unknown.
    If unknown, it will raise an error

    :param filename: path to file

    :return: current file type
    """
    with read_compressed_or_not(filename) as f:
        first_line = f.readline()
    if first_line.startswith("LOCUS       "):  # then this is probably a gbff/gbk file
        return "gbff"
    elif re.match(
        r"##gff-version\s{1,3}3", first_line
    ):  # prodigal gff header has two spaces between gff-version and 3... some gff user can have a tab
        return "gff"
    elif first_line.startswith(">"):
        return "fasta"
    elif "\t" in first_line:
        return "tsv"
    else:
        raise Exception(
            f"Filetype {filename} was not gff3 (file starts with '##gff-version 3') "
            "nor gbff/gbk (file starts with 'LOCUS       ') "
            "nor fasta (file starts with '>') nor tsv (file has '\t' in the first line). "
        )


def restricted_float(x: Union[int, float]) -> float:
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


def _plain_bfs(g: nx.Graph, source: Any, removed: set, weight: float):
    """
    A fast BFS node generator, copied from networkx then adapted to the current use case

    :param g: graph with the nodes
    :param source: current node
    :param removed: set of removed nodes
    :param weight: threshold to remove node or not

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
                        if (
                            len(edge_genes_n) / len(g.nodes[n]["genes"]) >= weight
                            and len(edge_genes_v) / len(g.nodes[v]["genes"]) >= weight
                        ):
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
        raise Exception(
            "If you give --clusters option, you must give at least --fasta or --anno"
        )

    if not any([args.fasta, args.anno]):
        raise Exception("At least one of --fasta or --anno must be given")

    if args.infer_singletons and args.clusters is None:
        logging.getLogger("PPanGGOLiN").warning(
            "--infer_singleton works only with --clusters."
        )


def parse_config_file(yaml_config_file: str) -> dict:
    """
    Parse yaml config file.

    :param yaml_config_file: config file in yaml

    :return: dict of config with key the command and as value another dict with param as key and value as value.
    """

    with yaml_config_file as yaml_fh:
        config = yaml.safe_load(yaml_fh)

    if config is None:
        config = {}

    # if config has a Parameters key. Update config with its content
    if config and "Parameters" in config:
        config.update(config["Parameters"])
        del config["Parameters"]

    # remove empty section that have no parameter specified in it. In this case they have a None value
    config = {
        section: param_val_dict
        for section, param_val_dict in config.items()
        if param_val_dict is not None
    }
    return config


def add_common_arguments(subparser: argparse.ArgumentParser):
    """
    Add common argument to the input subparser.

    :param subparser: A subparser object from any subcommand.
    """

    common = subparser._action_groups.pop(
        1
    )  # get the 'optional arguments' action group.
    common.title = "Common arguments"
    common.add_argument(
        "--verbose",
        required=False,
        type=int,
        default=1,
        choices=[0, 1, 2],
        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)",
    )
    common.add_argument(
        "--log",
        required=False,
        type=check_log,
        default="stdout",
        help="log output file",
    )
    common.add_argument(
        "-d",
        "--disable_prog_bar",
        required=False,
        action="store_true",
        help="disables the progress bars",
    )
    common.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force writing in output directory and in pangenome output file.",
    )
    common.add_argument(
        "--config",
        required=False,
        type=argparse.FileType(),
        help="Specify command arguments through a YAML configuration file.",
    )
    subparser._action_groups.append(common)


def get_arg_name(arg_val: Union[str, TextIOWrapper]) -> Union[str, TextIOWrapper]:
    """
    Returns the name of a file if the argument is a TextIOWrapper object,
    otherwise returns the argument value.

    :param arg_val: Either a string or a TextIOWrapper object.
    :return: Either a string or a TextIOWrapper object, depending on the type of the input argument.
    """

    if isinstance(arg_val, TextIOWrapper):
        return arg_val.name
    return arg_val


def overwrite_args(
    default_args: argparse.Namespace,
    config_args: argparse.Namespace,
    cli_args: argparse.Namespace,
):
    """
    Overwrite args objects.

    When arguments are given in CLI, their values are used instead of the ones found in the config file.
    When arguments are specified in the config file, they overwrite default values.

    :param default_args: default arguments
    :param config_args: arguments parsed from the config file
    :param cli_args: arguments parsed from the command line

    :return: final arguments
    """
    args = argparse.Namespace()
    all_params = [arg for arg in dir(default_args) if not arg.startswith("_")]

    for param in all_params:
        default_val = getattr(default_args, param)
        cli_val = getattr(cli_args, param, "unspecified")
        config_val = getattr(config_args, param, "unspecified")

        if param in cli_args and param not in config_args:
            # Use the value from the command line argument
            setattr(args, param, cli_val)

            if default_val != cli_val and param != "config":
                logging.getLogger("PPanGGOLiN").debug(
                    f'The parameter "--{param}: {get_arg_name(cli_val)}" has been specified in the command line with a non-default value.'
                    f" Its value overwrites the default value ({get_arg_name(default_val)})."
                )

        elif param not in cli_args and param in config_args:
            # Use the value from the config file
            setattr(args, param, config_val)

            if default_val != config_val:
                logging.getLogger("PPanGGOLiN").debug(
                    f'The parameter "--{param}: {get_arg_name(config_val)}" has been specified in the config file with a non-default value.'
                    f" Its value overwrites the default value ({get_arg_name(default_val)})."
                )

        elif param in cli_args and param in config_args:
            # Use the value from the command line argument (cli) if it's different from the config file (config)
            setattr(args, param, cli_val)

            if cli_val == config_val and cli_val != default_val:
                logging.getLogger("PPanGGOLiN").debug(
                    f'The parameter "--{param} {get_arg_name(cli_val)}" has been specified in both the command line '
                    f"and the config file with the same values, but with non-default value. "
                    f"Its value overwrites the default value ({get_arg_name(default_val)})."
                )

            elif cli_val != config_val and param != "config":
                # Values in cli and config differ. Use the value from the command line argument (cli)
                logging.getLogger("PPanGGOLiN").debug(
                    f'The parameter "--{param}" has been specified in both the command line ("{get_arg_name(cli_val)}") '
                    f'and the config file ("{get_arg_name(config_val)}") with different values. '
                    f"The value from the command line argument is used."
                )
        else:
            # Parameter is not defined in cli and in config. Use the default value.
            setattr(args, param, default_val)

    return args


def combine_args(args: argparse.Namespace, another_args: argparse.Namespace):
    """
    Combine two args object.

    :param args: initial arguments.
    :param another_args: another args

    :return: object with combined arguments
    """

    other_arg_names = [arg for arg in dir(another_args) if not arg.startswith("_")]

    for arg_name in other_arg_names:
        arg_val = getattr(another_args, arg_name)
        setattr(args, arg_name, arg_val)

    return args


def get_args_differing_from_default(
    default_args: argparse.Namespace,
    final_args: argparse.Namespace,
    param_to_ignore: Union[List[str], Set[str]] = None,
) -> dict:
    """
    Get the parameters that have different value than default values.

    :params default_args: default arguments
    :params final_args: final arguments to compare with default
    :params param_to_ignore: list of params to ignore.

    :return: A dict with param that differ from default as key and the final value of the param as value
    """
    param_to_ignore = [] if param_to_ignore is None else param_to_ignore
    all_params = [
        arg
        for arg in dir(final_args)
        if not arg.startswith("_")
        if arg not in param_to_ignore
    ]

    params_that_differ = {
        param: getattr(final_args, param)
        for param in all_params
        if getattr(default_args, param) != getattr(final_args, param)
    }

    return params_that_differ


def manage_cli_and_config_args(
    subcommand: str, config_file: str, subcommand_to_subparser: dict
) -> argparse.Namespace:
    """
    Manage command line and config arguments for the given subcommand.

    This function parse arguments from the cmd line and config file and set up the following priority: cli > config > default
    When the subcommand is a workflow, the subcommand used in workflows are also parsed in the config.


    :params subcommand: Name of the subcommand.
    :params config_file: Path to the config file given in argument. If None, only default and cli arguments value are used.
    :params subcommand_to_subparser: Dict with subcommand name as key and the corresponding subparser function as value.
    """
    if config_file:
        config = parse_config_file(config_file)
    else:
        config = {}

    # convert config dict to defaultdict
    config = defaultdict(dict, config)

    cmd_subparser = subcommand_to_subparser[subcommand]

    default_args = get_default_args(subcommand, cmd_subparser)

    cli_args = get_cli_args(cmd_subparser)

    all_cmd_param_names = {
        arg_name for arg_name in dir(default_args) if not arg_name.startswith("_")
    }

    input_params = {param for param in all_cmd_param_names if param in ALL_INPUT_PARAMS}

    general_params = {
        param for param in all_cmd_param_names if param in ALL_GENERAL_PARAMS
    }

    specific_params = all_cmd_param_names - (input_params | general_params)

    all_unspecific_params = ALL_INPUT_PARAMS + ALL_GENERAL_PARAMS
    # manage logging first to correctly set it up and to be able to log any issue when using config file later on
    config_general_args = get_config_args(
        subcommand,
        cmd_subparser,
        config,
        "general_parameters",
        general_params,
        strict_config_check=False,
    )
    general_args = overwrite_args(default_args, config_general_args, cli_args)

    set_verbosity_level(general_args)

    config_input_args = get_config_args(
        subcommand,
        cmd_subparser,
        config,
        "input_parameters",
        input_params,
        strict_config_check=True,
    )

    if subcommand in WORKFLOW_SUBCOMMANDS:
        # for workflow commands there is no section dedicated in the config: so no specific_args
        # only general_parameters and  sections of commands launched in the worklow commands are used
        config_args = combine_args(config_general_args, config_input_args)
    else:
        config_specific_args = get_config_args(
            subcommand,
            cmd_subparser,
            config,
            subcommand,
            specific_params,
            strict_config_check=True,
        )
        config_args = combine_args(config_general_args, config_specific_args)
        config_args = combine_args(config_args, config_input_args)

    # manage priority between source of args
    # cli > config > default

    args = overwrite_args(default_args, config_args, cli_args)
    params_that_differ = get_args_differing_from_default(
        default_args, args, input_params
    )

    if params_that_differ:
        params_that_differ_str = ", ".join(
            f"{p}={v}" for p, v in params_that_differ.items()
        )
        logging.getLogger("PPanGGOLiN").debug(
            f"{len(params_that_differ)} {subcommand} parameters have non-default value: {params_that_differ_str}"
        )

    # manage workflow command
    workflow_steps = []
    if subcommand in WORKFLOW_SUBCOMMANDS:

        for workflow_step in ALL_WORKFLOW_DEPENDENCIES:
            if (
                workflow_step in ["rgp", "spot"]
                and subcommand in ["workflow", "panmodule"]
            ) or (workflow_step == "module" and subcommand in ["workflow", "panrgp"]):
                continue
            logging.getLogger("PPanGGOLiN").debug(
                f"Parsing {workflow_step} arguments in config file."
            )
            step_subparser = subcommand_to_subparser[workflow_step]

            default_step_args = get_default_args(
                workflow_step, step_subparser, unwanted_args=all_unspecific_params
            )

            # remove general args
            all_param_names = {
                arg_name
                for arg_name in dir(default_step_args)
                if not arg_name.startswith("_")
            }
            specific_step_params = {
                param_name
                for param_name in all_param_names
                if param_name not in all_unspecific_params
            }
            config_step_args = get_config_args(
                workflow_step,
                step_subparser,
                config,
                workflow_step,
                specific_step_params,
                strict_config_check=True,
            )

            # overwrite write and draw default when not specified in config
            if workflow_step == "write_pangenome":
                for out_flag in WRITE_PAN_FLAG_DEFAULT_IN_WF:
                    if out_flag not in config[workflow_step]:
                        setattr(default_step_args, out_flag, True)

            if workflow_step == "write_genomes":
                for out_flag in WRITE_GENOME_FLAG_DEFAULT_IN_WF:
                    if out_flag not in config[workflow_step]:
                        setattr(default_step_args, out_flag, True)

            if workflow_step == "draw":
                for out_flag in DRAW_FLAG_DEFAULT_IN_WF:
                    if out_flag not in config[workflow_step]:
                        setattr(default_step_args, out_flag, True)

            step_args = overwrite_args(default_step_args, config_step_args, cli_args)

            step_params_that_differ = get_args_differing_from_default(
                default_step_args, step_args
            )

            if step_params_that_differ:
                step_params_that_differ_str = ", ".join(
                    f"{p}={v}" for p, v in step_params_that_differ.items()
                )
                logging.getLogger("PPanGGOLiN").debug(
                    f"{len(step_params_that_differ)} {workflow_step} parameters have "
                    f"a non-default value: {step_params_that_differ_str}"
                )

            # add step name to differentiate the params
            step_params_that_differ = {
                f"{workflow_step}:{param}": value
                for param, value in step_params_that_differ.items()
            }

            params_that_differ.update(step_params_that_differ)

            # Add args namespace of the step to the initial args namespace
            setattr(args, workflow_step, step_args)

    if params_that_differ:
        logging.getLogger("PPanGGOLiN").info(
            f"{len(params_that_differ)} parameters have a non-default value."
        )

    check_config_consistency(config, workflow_steps)

    return args


def check_config_consistency(config: dict, workflow_steps: list):
    """
    Check that the same parameter used in different subcommand inside a workflow has the same value.

    If not, the function throw a logging.getLogger("PPanGGOLiN").warning.

    :params config_dict: config dict with as key the section of the config file and as value another dict pairing name and value of parameters.
    :params workflow_steps: list of subcommand names used in the workflow execution.
    """

    def count_different_values(values: Iterable[Union[int, str, Tuple, List]]) -> int:
        """
        Returns the number of unique values in a list.

        :param values: A list of values to count.
        :return: The number of unique values in the list.
        """
        hashable_values = set()
        for value in values:
            hashable_value = tuple(value) if isinstance(value, list) else value
            hashable_values.add(hashable_value)
        return len(hashable_values)

    # params used in multiple subcommands
    all_params = [
        param
        for subcmd, param_to_value_dict in config.items()
        for param in param_to_value_dict
        if subcmd in workflow_steps
    ]
    duplicate_params = [param for param in all_params if all_params.count(param) > 1]

    for duplicate_param in set(duplicate_params):
        step_to_value = {
            step: param_to_value[duplicate_param]
            for step, param_to_value in config.items()
            if duplicate_param in param_to_value
        }

        if count_different_values(step_to_value.values()) > 1:
            logging.getLogger("PPanGGOLiN").warning(
                f"The parameter {duplicate_param} used in multiple subcommands of the workflow is specified with different values in config file: {step_to_value}."
            )


def set_up_config_param_to_parser(config_param_val: dict) -> list:
    """
    Take dict pairing parameters and values and format the corresponding list of arguments to feed a parser.

    When the parameter value is False, the parameter is a flag and thus is not added to the list.

    :params config_param_val: Dict with parameter name as key and parameter value as value.

    :return: list of argument strings formatted for an argparse.ArgumentParser object.
    """

    arguments_to_parse = []
    for param, val in config_param_val.items():

        if isinstance(val, bool) or val is None or val == "None":
            # param is a flag
            if val is True:
                arguments_to_parse.append(f"--{param}")
            # if val is False or None we don't add id to the
        else:
            arguments_to_parse.append(f"--{param}")

            if isinstance(val, list):
                # range of values need to be added one by one
                arguments_to_parse += [str(v) for v in val]
            else:
                arguments_to_parse.append(str(val))
    return arguments_to_parse


def get_subcommand_parser(
    subparser_fct: Callable, name: str = ""
) -> Tuple[argparse._SubParsersAction, argparse.ArgumentParser]:
    """
    Get subcommand parser object using the given subparser function.

    Common arguments are also added to the parser object.

    :params subparser_fct:
    :name: Name of section to add more info in the parser in case of error.

    :return: The parser and subparser objects
    """
    prog = ""
    usage = ""

    if name:
        prog = f"Parsing section {name} in config file"
        usage = "Yaml config file"

    parser = argparse.ArgumentParser(prog=prog, allow_abbrev=False, add_help=False)

    subparsers = parser.add_subparsers(
        metavar="", dest="subcommand", title="subcommands", description=""
    )

    sub = subparser_fct(subparsers)
    sub.usage = usage
    add_common_arguments(sub)

    # set off required flag in required arguments
    for arg_action in sub._actions:
        if arg_action.required:
            arg_action.required = False
    return parser, sub


def get_default_args(
    subcommand: str, subparser_fct: Callable, unwanted_args: list = None
) -> argparse.Namespace:
    """
    Get default value for the arguments for the given subparser function.

    :params subcommand: Name of the ppanggolin subcommand.
    :params subparser_fct: Subparser function to use. This subparser give the expected argument for the subcommand.
    :params unwanted_args: List of arguments to filter out.

    :return args: arguments with default values.
    """
    unwanted_args = [] if unwanted_args is None else unwanted_args
    parser, sub = get_subcommand_parser(subparser_fct, subcommand)

    # remove unwanted argumnents
    sub._actions = [
        p_action for p_action in sub._actions if p_action.dest not in unwanted_args
    ]

    args = parser.parse_args([subcommand])

    return args


def get_config_args(
    subcommand: str,
    subparser_fct: Callable,
    config_dict: dict,
    config_section: str,
    expected_params: Union[List[str], Set[str]],
    strict_config_check: bool,
) -> argparse.Namespace:
    """
    Parsing parameters of a specific section of the config file.

    If some parameter are not specified in the config they are not added to the args object.

    :params subcommand: Name of the ppanggolin subcommand.
    :params subparser_fct: Subparser function to use. This subparser give the expected argument for the subcommand.
    :params config_dict: config dict with as key the section of the config file and as value another dict pairing name and value of parameters.
    :params config_section: Which section to parse in config file.
    :params expected_params: List of argument to expect in the parser. If the parser has other arguments, these arguments are filtered out.
    :params strict_config_check: if set to true, an error is raised when a parameter is found in the config which it is not in the expected_params list.

    :return args: Arguments parse from the config
    """
    config = config_dict[config_section]

    parser, sub = get_subcommand_parser(subparser_fct, subcommand)

    # for all args set default to None to be able to distinguish params that have been specified in config
    erase_default_value(sub)

    # Manage args
    sub._actions = [
        p_action for p_action in sub._actions if p_action.dest in expected_params
    ]

    if not strict_config_check:
        # remove param found in config that are not expected by parser. useful for general_parameters.
        expected_args_names = [p_action.dest for p_action in sub._actions]
        unexpected_config = [
            f"{name}:{value}"
            for name, value in config.items()
            if name not in expected_args_names
        ]
        config = {
            name: value for name, value in config.items() if name in expected_args_names
        }

        if unexpected_config:
            logging.getLogger("PPanGGOLiN").info(
                f"While parsing {config_section} section in config file, {len(unexpected_config)} unexpected parameters "
                f'were ignored : {" ".join(unexpected_config)}'
            )
    else:
        for param_name in config:
            if param_name not in expected_params:
                sub.error(f"unrecognized arguments: {param_name}")

    config_args_to_parse = set_up_config_param_to_parser(config)

    args = parser.parse_args([subcommand] + config_args_to_parse)

    # remove argument that have not been specified in config file
    # unspecified argument have None as value
    delete_unspecified_args(args)

    return args


def get_cli_args(subparser_fct: Callable) -> argparse.Namespace:
    """
    Parse command line arguments using the specified parsing function.

    :params subparser_fct: Subparser function to use. This subparser give the expected argument for the subcommand.
    """

    parser, sub = get_subcommand_parser(subparser_fct)

    # for all args set default to None to be able to distinguish params that have been specified in config
    erase_default_value(sub)

    cli_args = parser.parse_args()  # parse cli

    # remove argument that have not been specified
    delete_unspecified_args(cli_args)
    delattr(cli_args, "subcommand")

    return cli_args


def erase_default_value(parser: argparse.ArgumentParser):
    """
    Remove default action in the given list of argument parser actions.

    This is dnoe to distinguish specified arguments.

    :params parser: An argparse.ArgumentParser object with default values to erase.
    """

    # for all args set default to None
    for p_action in parser._actions:
        p_action.default = None


def delete_unspecified_args(args: argparse.Namespace):
    """
    Delete argument from the given argparse.Namespace with None values.

    :param args: arguments to filter.
    """

    for arg_name, arg_val in args._get_kwargs():
        if arg_val is None:
            delattr(args, arg_name)


def extract_contig_window(
    contig_size: int,
    positions_of_interest: Iterable[int],
    window_size: int,
    is_circular: bool = False,
):
    """
    Extracts contiguous windows around positions of interest within a contig.

    :param contig_size: Number of genes in contig.
    :param positions_of_interest: An iterable containing the positions of interest.
    :param window_size: The size of the window to extract around each position of interest.
    :param is_circular: Indicates if the contig is circular.
    :return: Yields tuples representing the start and end positions of each contiguous window.
    """
    windows_coordinates = []

    # Sort the positions of interest
    sorted_positions = sorted(positions_of_interest)

    # Check if any position of interest is out of range
    if sorted_positions[0] < 0 or sorted_positions[-1] >= contig_size:
        raise IndexError(
            f"Positions of interest are out of range. "
            f"Contig has {contig_size} genes while given min={sorted_positions[0]} & max={sorted_positions[-1]} positions"
        )

    if is_circular:
        first_position = sorted_positions[0]
        last_position = sorted_positions[-1]
        # in a circular contig, if the window of a gene of interest overlaps the end/start of the contig
        # an out of scope position is added to the sorted positions to take into account those positions
        # the returned window are always checked that its positions are not out of range...
        # so there's no chance to find an out of scope position in final list
        if first_position - window_size < 0:
            out_of_scope_position = contig_size + first_position
            sorted_positions.append(out_of_scope_position)

        if last_position + window_size >= contig_size:
            out_of_scope_position = last_position - contig_size
            sorted_positions.insert(0, out_of_scope_position)

    start_po = max(sorted_positions[0] - window_size, 0)

    for position, next_po in zip_longest(sorted_positions, sorted_positions[1:]):

        if next_po is None:
            # If there are no more positions, add the final window
            end_po = min(position + window_size, contig_size - 1)
            windows_coordinates.append((start_po, end_po))

        elif position + window_size + 1 < next_po - window_size:
            # If there is a gap between positions, add the current window
            # and update the start position for the next window
            end_po = min(position + window_size, contig_size - 1)

            windows_coordinates.append((start_po, end_po))

            start_po = max(next_po - window_size, 0)

    return windows_coordinates


def parse_input_paths_file(
    path_list_file: Path,
) -> Dict[str, Dict[str, Union[Path, List[str]]]]:
    """
    Parse an input paths file to extract genome information.

    This function reads an input paths file, which is in TSV format, and extracts genome information
    including file paths and putative circular contigs.

    :param path_list_file: The path to the input paths file.
    :return: A dictionary where keys are genome names and values are dictionaries containing path information and
             putative circular contigs.
    :raises FileNotFoundError: If a specified genome file path does not exist.
    :raises Exception: If there are no genomes in the provided file.
    """
    logging.getLogger("PPanGGOLiN").info(
        f"Reading {path_list_file} to process genome files"
    )
    genome_name_to_genome_path = {}

    for line in read_compressed_or_not(path_list_file):
        elements = [el.strip() for el in line.split("\t")]
        genome_file_path = Path(elements[1])
        genome_name = elements[0]
        putative_circular_contigs = elements[2:]

        if not genome_file_path.exists():
            # Check if the file path doesn't exist and try an alternative path.
            genome_file_path_alt = path_list_file.parent.joinpath(genome_file_path)

            if not genome_file_path_alt.exists():
                raise FileNotFoundError(
                    f"The file path '{genome_file_path}' for genome '{genome_name}' specified in '{path_list_file}' does not exist."
                )
            else:
                genome_file_path = genome_file_path_alt

        genome_name_to_genome_path[genome_name] = {
            "path": genome_file_path,
            "circular_contigs": putative_circular_contigs,
        }

    if len(genome_name_to_genome_path) == 0:
        raise Exception(f"There are no genomes in the provided file: {path_list_file} ")

    return genome_name_to_genome_path


def flatten_nested_dict(
    nested_dict: Dict[str, Union[Dict, int, str, float]],
) -> Dict[str, Union[int, str, float]]:
    """
    Flattens a nested dictionary into a flat dictionary by concatenating keys at different levels.

    :param nested_dict: The nested dictionary to be flattened.
    :return: A flat dictionary with concatenated keys.
    """
    flat_dict = {}

    def flatten(dictionary, parent_key=""):
        for key, val in dictionary.items():
            new_key = f"{parent_key}_{key}" if parent_key else key
            if isinstance(val, dict):
                flatten(val, new_key)
            else:
                flat_dict[new_key] = val

    flatten(nested_dict)
    return flat_dict


def get_major_version(version: str) -> int:
    """
    Extracts the major version number from a version string.

    :param version: A string representing the version number.
    :return: The major version extracted from the input version string.
    :raises ValueError: If the input version does not have the expected format.
    """
    try:
        major_version = int(version.split(".")[0])
    except ValueError:
        raise ValueError(f"Version {version} does not have the expected format.")

    return major_version


def check_version_compatibility(file_version: str) -> None:
    """
    Checks the compatibility of the provided pangenome file version with the current PPanGGOLiN version.

    :param file_version: A string representing the version of the pangenome file.
    """
    # Get the current PPanGGOLiN version
    current_version = distribution("ppanggolin").version

    current_version_major = get_major_version(current_version)
    file_major_version = get_major_version(file_version)

    # Check for compatibility issues
    if file_major_version != current_version_major:
        logging.getLogger("PPanGGOLiN").error(
            "Your pangenome file has been created with a different major version "
            "of PPanGGOLiN than the one installed in the system. This mismatch may lead to compatibility issues."
        )

    if file_major_version < 2 and current_version_major >= 2:
        raise ValueError(
            f"The provided pangenome file was created by PPanGGOLiN version {file_version}, which is "
            f"incompatible with the current PPanGGOLiN version {current_version}."
        )


def find_consecutive_sequences(sequence: List[int]) -> List[List[int]]:
    """
    Find consecutive sequences in a list of integers.

    :param sequence: The input list of integers.

    :return: A list of lists containing consecutive sequences of integers.
    """
    s_sequence = sorted(sequence)

    consecutive_sequences = [[s_sequence[0]]]

    for index in s_sequence[1:]:
        if index == consecutive_sequences[-1][-1] + 1:
            consecutive_sequences[-1].append(index)
        else:
            # there is a break in the consecutivity
            consecutive_sequences.append([index])

    return consecutive_sequences


def find_region_border_position(
    region_positions: List[int], contig_gene_count: int
) -> Tuple[int, int]:
    """
    Find the start and stop integers of the region considering circularity of the contig.

    :param region_positions: List of positions that belong to the region.
    :param contig_gene_count: Number of gene in the contig. The contig is considered circular.

    :return: A tuple containing the start and stop integers of the region.
    """

    consecutive_region_positions = get_consecutive_region_positions(
        region_positions, contig_gene_count
    )

    return consecutive_region_positions[0][0], consecutive_region_positions[-1][-1]


def get_consecutive_region_positions(
    region_positions: List[int], contig_gene_count: int
) -> List[List[int]]:
    """
    Order integers position of the region considering circularity of the contig.

    :param region_positions: List of positions that belong to the region.
    :param contig_gene_count: Number of gene in the contig. The contig is considered circular.

    :return: An ordered list of integers of the region.

    :raises ValueError: If unexpected conditions are encountered.
    """
    if len(region_positions) == 0:
        raise ValueError("Region has no position. This is unexpected.")

    consecutive_sequences = sorted(find_consecutive_sequences(region_positions))

    if len(consecutive_sequences) == 0:
        raise ValueError(
            "No consecutive sequences found in the region. This is unexpected."
        )

    elif len(consecutive_sequences) == 1:
        return consecutive_sequences

    elif len(consecutive_sequences) == 2:
        # Check for overlaps at the edge of the contig
        if consecutive_sequences[0][0] != 0:
            raise ValueError(
                f"Two sequences of consecutive positions ({consecutive_sequences}) "
                f"indicate an overlap on the edge of the contig, but neither starts at the beginning of the contig (0)."
            )

        elif consecutive_sequences[-1][-1] != contig_gene_count - 1:
            raise ValueError(
                f"Two sequences of consecutive positions ({consecutive_sequences}) "
                f"indicate an overlap on the edge of the contig, but neither ends at the end of the contig ({contig_gene_count - 1})."
            )

        return [consecutive_sequences[-1], consecutive_sequences[0]]

    elif len(consecutive_sequences) > 2:
        raise ValueError(
            f"More than two consecutive sequences found ({len(consecutive_sequences)}). "
            f"This is unexpected. Consecutive sequences: {consecutive_sequences}. "
            "The region should consist of consecutive positions along the contig."
        )


def run_subprocess(
    cmd: List[str],
    output: Optional[Path] = None,
    msg: str = "Subprocess failed with the following error:\n",
):
    """
    Run a subprocess command and write the output to the given path.

    :param cmd: List of program arguments.
    :param output: Path to write the subprocess output (optional).
    :param msg: Message to print if the subprocess fails.

    :raises FileNotFoundError: If the command's executable is not found.
    :raises subprocess.CalledProcessError: If the subprocess returns a non-zero exit code.
    """
    if not cmd or shutil.which(cmd[0]) is None:
        raise FileNotFoundError(
            f"Command '{cmd[0]}' not found. Please install it and try again."
        )

    logging.getLogger("PPanGGOLiN").debug(f"Running command: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as subprocess_err:
        logging.getLogger("PPanGGOLiN").error(msg)
        for line in subprocess_err.stderr.strip().split("\n"):
            logging.getLogger("PPanGGOLiN").error(line)
        raise
    else:
        if output:
            with open(output, "w") as fout:
                fout.write(result.stdout)


def has_non_ascii(string_to_test: Union[str, Collection[str]]) -> bool:
    """
    Check if a string or a collection of strings contains any non-ASCII characters.

    :param string_to_test: A single string or a collection of strings to check.
    :return: True if any string contains non-ASCII characters, False otherwise.
    """
    try:
        if isinstance(string_to_test, str):
            string_to_test.encode("ascii")
        else:
            for item in string_to_test:
                item.encode("ascii")
        return False
    except UnicodeEncodeError:
        return True


def replace_non_ascii(
    string_with_ascii: Union[str, Collection[str]], replacement_string: str = "_"
) -> Union[str, Collection[str]]:
    """
    Replace all non-ASCII characters in a string or a collection of strings
    with a specified replacement string.

    :param string_with_ascii: A string or collection of strings potentially containing non-ASCII characters.
    :param replacement_string: The string to replace non-ASCII characters with (default is '_').
    :return: A new string or collection where all non-ASCII characters have been replaced.
    """

    def replace(s: str) -> str:
        return re.sub(r"[^\x00-\x7F]+", replacement_string, s)

    if isinstance(string_with_ascii, str):
        return replace(string_with_ascii)
    return type(string_with_ascii)(replace(s) for s in string_with_ascii)


def check_tools_availability(
    tool_to_description: Union[Dict[str, str], List[str]],
) -> dict[str, bool]:
    """
    Check if the given command-line tools are available in the system's PATH.

    :param tool_to_description: A dictionary where keys are tool names and values are descriptions of their purpose, or a list of tool names.
    :return: A dictionary with tool names as keys and boolean values indicating availability.
    """
    if isinstance(tool_to_description, list):
        tool_to_description = {tool: "" for tool in tool_to_description}

    availability = {
        tool: shutil.which(tool) is not None for tool in tool_to_description
    }

    for tool, is_available in availability.items():
        if not is_available:
            caller_frame = inspect.stack()[1]
            caller_module = caller_frame.frame.f_globals["__name__"]
            caller_function = caller_frame.function

            logging.getLogger("PPanGGOLiN").warning(
                f"Missing required command: '{tool}' {tool_to_description[tool]}. "
                f"This check was triggered in '{caller_function}' inside module '{caller_module}'. "
                "Please install the missing command to ensure proper functionality of PPanGGOLiN."
            )

    return availability
