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
from ppanggolin.annotate import annotatePangenome, readAnnotations, getGeneSequencesFromFastas
from ppanggolin.cluster import clustering, readClustering
from ppanggolin.graph import computeNeighborsGraph
from ppanggolin.evolution import makeEvolutionCurve
from ppanggolin.partition import partition
from ppanggolin.formats import writePangenome, readPangenome
### a global workflow that does everything in one go.


def launch(args):

    pangenome = Pangenome()
    logging.getLogger().info(f"Starting RAM : {getCurrentRAM()}")
    filename = mkFilename(args.basename, args.output, args.force)
    if args.gff:#if the gff is provided, we read annotations from it
        getSeq = True
        if args.clusters is not None:
            getSeq = False
        readAnnotations(pangenome, args.gff, getSeq)
        logging.getLogger().info(f"post reading annotations : {getCurrentRAM()}")
        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is None:
            raise Exception("The gff provided did not have any sequence informations, you did not provide clusters and you did not provide fasta file. Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is not None:
            getGeneSequencesFromFastas(pangenome, args.fasta)
        
        if args.clusters is not None:
            readClustering(pangenome, args.clusters)

        elif args.clusters is None:#we should have the sequences here.
            clustering(pangenome, args.tmpdir, args.cpu)
        logging.getLogger().info(f"post clustering : {getCurrentRAM()}")
    elif args.fasta is not None:
            logging.getLogger().info("You did not provide annotations but provided fasta. Everything will be done from here using the fasta only.")
            pangenome = Pangenome()
            annotatePangenome(pangenome, args.fasta, args.tmpdir, args.cpu)
            clustering(pangenome, args.tmpdir, args.cpu)
    computeNeighborsGraph(pangenome)
    logging.getLogger().info(f"post making the graph : {getCurrentRAM()}")
    partition(pangenome, cpu = args.cpu, tmpdir = args.tmpdir)
    logging.getLogger().info(f"post partitionning : {getCurrentRAM()}")
    makeEvolutionCurve(pangenome,args.output, args.tmpdir, cpu=args.cpu)
    writePangenome(pangenome, filename, args.force)



def workflowSubparser(subparser):
    parser = subparser.add_parser("workflow",help = "Easy workflow to run a pangenome analysis in one go without parameter tuning")
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--basename",required = False, default = "pangenome", help = "basename for the output file")
    required = parser.add_argument_group(title = "Input arguments", description = "The possible input arguments :")
    required.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism. This option can be used alone.")
    required.add_argument('--gff', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. This option can be used alone IF the fasta files are in the gff files, otherwise --fasta needs to be used.")
    required.add_argument("--clusters",required=False, type=str, help = "a tab-separated file listing the cluster names, the gene IDs, and optionnally whether they are a fragment or not.")
    return parser