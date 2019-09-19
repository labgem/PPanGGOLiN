#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import os
import time
import argparse

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mkFilename
from ppanggolin.annotate import annotatePangenome, readAnnotations, getGeneSequencesFromFastas
from ppanggolin.cluster import clustering, readClustering
from ppanggolin.graph import computeNeighborsGraph
from ppanggolin.nem.evolution import makeEvolutionCurve
from ppanggolin.nem.partition import partition
from ppanggolin.formats import writePangenome, writeFlatFiles
from ppanggolin.figures import drawTilePlot, drawUCurve
### a global workflow that does everything in one go.

def launch(args):
    pangenome = Pangenome()
    filename = mkFilename(args.basename, args.output, args.force)
    if args.anno:#if the annotations are provided, we read from it
        getSeq = True
        if args.clusters is not None:
            getSeq = False
        readAnnotations(pangenome, args.anno, getSeq)
        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is None:
            raise Exception("The gff/gbff provided did not have any sequence informations, you did not provide clusters and you did not provide fasta file. Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is not None:
            getGeneSequencesFromFastas(pangenome, args.fasta)

        if args.clusters is not None:
            readClustering(pangenome, args.clusters)

        elif args.clusters is None:#we should have the sequences here.
            clustering(pangenome, args.tmpdir, args.cpu)
    elif args.fasta is not None:
            pangenome = Pangenome()
            annotatePangenome(pangenome, args.fasta, args.tmpdir, args.cpu)
            clustering(pangenome, args.tmpdir, args.cpu)

    computeNeighborsGraph(pangenome)
    partition(pangenome, tmpdir = args.tmpdir, cpu = args.cpu)
    if args.evol:
        makeEvolutionCurve(pangenome,args.output, args.tmpdir, cpu=args.cpu)

    drawTilePlot(pangenome, args.output, nocloud = False if len(pangenome.organisms) < 500 else True )
    drawUCurve(pangenome, args.output)

    writeFlatFiles(pangenome, args.output, args.cpu, csv = True, genePA=True, gexf=True, light_gexf = True, projection=True, json = True, stats = True, partitions = True)

    writePangenome(pangenome, filename, args.force)

def workflowSubparser(subparser):
    parser = subparser.add_parser("workflow", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Input arguments", description = "The possible input arguments :")
    required.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism. This option can be used alone.")
    required.add_argument('--anno', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used.")
    required.add_argument("--clusters",required=False, type=str, help = "a tab-separated file listing the cluster names, the gene IDs, and optionnally whether they are a fragment or not.")

    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--basename",required = False, default = "pangenome", help = "basename for the output file")
    optional.add_argument("--evol", required=False, action = "store_true", help = "Use to compute the evolution curves (WARNING: can be time consumming)")

    return parser