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
from numpy import repeat


def check_log(name):
    if name == "stdout":
        return sys.stdout
    elif name == "stderr":
        return sys.stderr
    else:
        return open(name, "w")


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
        reads a file object or file path, uncompresses it if need be.
        returns a TextIO object in read only.
    """
    file = file_or_file_path
    if isinstance(file, str):
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except AttributeError:
            return (file)
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        return (TextIOWrapper(gzip.open(filename=file, mode="r")))
    else:
        file.close()
        file = open(file.name, "r")
        return (file)


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


def mkOutdir(output, force):
    if not os.path.exists(output):
        os.makedirs(output)
    elif not force:
        raise FileExistsError(f"{output} already exists. Use -f if you want to overwrite the files in the directory")


def mkFilename(basename, output, force):
    """
        Returns a usable filename for a ppanggolin output file, or crashes.
    """
    filename = Path(output + "/" + basename)
    if filename.suffix != ".h5":
        filename = filename.with_suffix(".h5")

    mkOutdir(output, force)

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
                obj["genes"][gene.family] = set([gene])
            except KeyError:
                obj["genes"] = {gene.family: set([gene])}
    else:
        try:
            obj["genes"].add(gene)
        except KeyError:
            obj["genes"] = set([gene])


def check_option_workflow(args):
    if args.clusters is not None and not any([args.fasta, args.anno]):
        raise Exception("If you give --clusters option, you must give at least --fasta or --anno")

    if not any([args.fasta, args.anno]):
        raise Exception("At least one of --fasta or --anno must be given")