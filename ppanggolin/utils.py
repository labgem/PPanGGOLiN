#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys
import os
import gzip
import mmap
import argparse
from io import TextIOWrapper
from pathlib import Path

import pkg_resources
from numpy import repeat
import logging


def check_log(name):
    if name == "stdout":
        return sys.stdout
    elif name == "stderr":
        return sys.stderr
    else:
        return open(name, "w")


def check_tsv_sanity(tsv):
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
        raise Exception(
            f"Some of your genomes have identical names. The duplicated names are the following : "
            f"'{' '.join(duplicated_names)}'")


def check_input_files(anno=None, pangenome=None, fasta=None):
    """
        Checks if the provided input files exist and are of the proper format
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


def jaccard_similarities(mat, jaccard_similarity_th):
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


def read_compressed_or_not(file_or_file_path):
    """
        reads a file object or file path, uncompresses it, if need be.
        returns a TextIO object in read only.
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


def write_compressed_or_not(file_path, compress):
    """
        Returns a file-like object, compressed or not.
    """
    if compress:
        return gzip.open(file_path + ".gz", mode="wt")
    else:
        return open(file_path, "w")


def is_compressed(file_or_file_path):
    """
        Checks is a file, or file path given is compressed or not
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


def get_num_lines(file):
    fp = open(file, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


def mk_outdir(output, force):
    if not os.path.exists(output):
        os.makedirs(output)
    elif not force:
        raise FileExistsError(f"{output} already exists. Use -f if you want to overwrite the files in the directory")


def mk_file_name(basename, output, force):
    """
        Returns a usable filename for a ppanggolin output file, or crashes.
    """
    filename = Path(output + "/" + basename)
    if filename.suffix != ".h5":
        filename = filename.with_suffix(".h5")

    mk_outdir(output, force)

    if filename.exists() and not force:
        raise FileExistsError(f"{filename.name} already exists. Use -f if you want to overwrite the file")
    return filename


def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def min_one(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("%r is inferior to 1" % (x,))
    return x


def connected_components(g, removed, weight):
    """
        Yields subgraphs of each connected component you get when filtering edges based on the given weight.
    """
    for v in g.nodes:
        if v not in removed:
            c = set(_plain_bfs(g, v, removed, weight))
            yield c
            removed.update(c)


def _plain_bfs(g, source, removed, weight):
    """A fast BFS node generator, copied from networkx then adapted to the current use case"""
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


def add_gene(obj, gene, fam_split=True):
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
    if args.clusters is not None and not any([args.fasta, args.anno]):
        raise Exception("If you give --clusters option, you must give at least --fasta or --anno")

    if not any([args.fasta, args.anno]):
        raise Exception("At least one of --fasta or --anno must be given")
