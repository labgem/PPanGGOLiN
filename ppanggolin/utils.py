#!/usr/bin/env python3
#coding:utf-8

#default libraries
import sys
import gzip
from io import TextIOWrapper
import mmap
from pathlib import Path
import logging
import os

#installed libraries
import psutil

def read_compressed_or_not(file_or_file_path):
    """
        reads a file object or file path, uncompresses it if need be.
        returns a TextIO object in read only.
    """
    file = file_or_file_path
    if type(file) == str:
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except:
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
    if type(file) == str:
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except:
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

def getCurrentRAM():
    units = ["o","Ko","Mo","Go","To"]
    mem = float(psutil.virtual_memory()._asdict()["used"])
    unit = 0
    while mem >= 1024:
        mem = mem / 1024
        unit +=1
    return str(round(mem,3)) + " " + units[unit]

def mkFilename(basename, output, force):
    """
        Returns a usable filename for a ppanggolin output file, or crashes.
    """
    filename = Path(output + "/" + basename )
    if filename.suffix != ".h5":
        filename = filename.with_suffix(".h5")
    
    if not os.path.exists(output):
        os.makedirs(output)
    elif filename.exists() and not force:
        logging.getLogger().error(f"{filename.name} already exists. Use -f if you want to overwrite the file")
        exit(1)
    return filename