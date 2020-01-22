#!/usr/bin/env python3
#coding:utf-8

#default libraries
import os
import time
import argparse

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mkFilename
from ppanggolin.annotate import annotatePangenome, readAnnotations, getGeneSequencesFromFastas
from ppanggolin.cluster import clustering, readClustering
from ppanggolin.graph import computeNeighborsGraph
from ppanggolin.nem.rarefaction import makeRarefactionCurve
from ppanggolin.nem.partition import partition
from ppanggolin.formats import writePangenome, writeFlatFiles
from ppanggolin.figures import drawTilePlot, drawUCurve
from ppanggolin.info import printInfo
### a global workflow that does everything in one go.

def launch(args):
    pangenome = Pangenome()
    filename = mkFilename(args.basename, args.output, args.force)
    if args.anno:#if the annotations are provided, we read from it
        getSeq = True
        if args.clusters is not None:
            getSeq = False
        readAnnotations(pangenome, args.anno, getSeq)
        writePangenome(pangenome, filename, args.force)
        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta_list is None and args.fastas is None:
            raise Exception("The gff/gbff provided did not have any sequence informations, you did not provide clusters and you did not provide fasta file. Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta_list is not None:
            getGeneSequencesFromFastas(pangenome, args.fasta_list)

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fastas is not None:
            getGeneSequencesFromFastas(pangenome, args.fastas)

        if args.clusters is not None:
            readClustering(pangenome, args.clusters)

        elif args.clusters is None:#we should have the sequences here.
            clustering(pangenome, tmpdir = args.tmpdir, cpu = args.cpu, defrag = args.defrag)
    elif args.fasta_list is not None:
        annotatePangenome(pangenome, args.fasta_list, args.tmpdir, args.cpu, fastaList=True)
        writePangenome(pangenome, filename, args.force)
        clustering(pangenome, tmpdir = args.tmpdir,cpu = args.cpu, defrag = args.defrag)
    elif args.fastas is not None:
        annotatePangenome(pangenome, args.fastas, args.tmpdir, args.cpu, fastaList=False)
        writePangenome(pangenome, filename, args.force)
        clustering(pangenome, tmpdir = args.tmpdir,cpu = args.cpu, defrag = args.defrag)

    computeNeighborsGraph(pangenome)

    partition(pangenome, tmpdir = args.tmpdir, cpu = args.cpu, K=args.nb_of_partitions)
    writePangenome(pangenome, filename, args.force)

    if args.rarefaction:
        makeRarefactionCurve(pangenome,args.output, args.tmpdir, cpu=args.cpu)
    if len(pangenome.organisms) > 1 and len(pangenome.organisms) < 5000:
        drawTilePlot(pangenome, args.output, nocloud = False if len(pangenome.organisms) < 500 else True)
    drawUCurve(pangenome, args.output)

    writeFlatFiles(pangenome, args.output, args.cpu, csv = True, genePA=True, gexf=True, light_gexf = True, projection=True, json = True, stats = True, partitions = True)

    printInfo(filename, content = True)

def workflowSubparser(subparser):
    parser = subparser.add_parser("workflow", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Input arguments", description = "The possible input arguments :")
    required.add_argument('--fasta_list',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism. This option can be used alone. If this option is used the option --fastas is ignored")
    required.add_argument('--fastas',  required=False, nargs='+', type=str, help="All the fasta files to be used, the filepaths will be used as organism names.")
    required.add_argument('--anno', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. This option can be used alone IF the fasta sequences are in the gff files, otherwise the argument --fasta_list or --fastas needs to be used.")
    required.add_argument("--clusters",required=False, type=str, help = "a tab-separated file listing the cluster names, the gene IDs, and optionnally whether they are a fragment or not.")

    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--basename",required = False, default = "pangenome", help = "basename for the output file")
    optional.add_argument("--rarefaction", required=False, action = "store_true", help = "Use to compute the rarefaction curves (WARNING: can be time consumming)")
    optional.add_argument("-K","--nb_of_partitions",required=False, default=-1, type=int, help = "Number of partitions to use. Must be at least 3. If under 3, it will be detected automatically.")
    optional.add_argument("--defrag",required=False, action="store_true", help = "Realign gene families to associated fragments with their non-fragmented gene family.")
    return parser
