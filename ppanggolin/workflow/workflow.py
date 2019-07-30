#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import logging
import os
import time

#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import getCurrentRAM, mkFilename
from ppanggolin.annotate import annotatePangenome
from ppanggolin.cluster import clustering
from ppanggolin.formats import writePangenome
### a global workflow that does everything in one go.


def launch(args):
    if args.fasta is not None:
        if args.gff is not None:
            #read annotations
            raise NotImplementedError()
        else:
            #do everything.
            filename = mkFilename(args.basename, args.output, args.force)
            pangenome = Pangenome()
            annotatePangenome(pangenome, args.fasta, args.tmpdir, args.cpu)
            clustering(pangenome, args.tmpdir, args.cpu)
            #makegraph
            #partition pangenome
            
            writePangenome(pangenome, filename, args.force)



def workflowSubparser(subparser):
    parser = subparser.add_parser("workflow",help = "Easy workflow to run a pangenome analysis in one go without parameter tuning")
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--basename",required = False, default = "pangenome", help = "basename for the output file")
    required = parser.add_argument_group(title = "Input arguments", description = "The possible input arguments :")
    required.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism. This option can be used alone.")
    required.add_argument('--gff', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. This option can be used alone IF the fasta files are in the gff files, otherwise --fasta needs to be used.")
    return parser