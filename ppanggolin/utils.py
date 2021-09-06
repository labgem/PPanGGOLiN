#!/usr/bin/env python3
#coding:utf-8

#default libraries
import gzip
from io import TextIOWrapper
import mmap
from pathlib import Path
import os
import numpy
import argparse

def jaccard_similarities(mat,jaccard_similarity_th):
    cols_sum = mat.getnnz(axis=0)
    ab = mat.T * mat
    # for rows
    aa = numpy.repeat(cols_sum, ab.getnnz(axis=0))
    # for columns
    bb = cols_sum[ab.indices]
    similarities = ab.copy()
    similarities.data /= (aa + bb - ab.data)
    similarities.data[similarities.data<jaccard_similarity_th] = 0
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
            return(file)
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        return(TextIOWrapper(gzip.open(filename=file, mode="r")))
    else:
        file.close()
        file = open(file.name, "r")
        return(file)

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
    filename = Path(output + "/" + basename )
    if filename.suffix != ".h5":
        filename = filename.with_suffix(".h5")

    mkOutdir(output, force)

    if filename.exists() and not force:
        raise FileExistsError(f"{filename.name} already exists. Use -f if you want to overwrite the file")
    return filename

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def min_one(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("%r is inferior to 1"%(x,))
    return x