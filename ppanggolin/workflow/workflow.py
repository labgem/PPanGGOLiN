#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import time
import argparse

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mkFilename, min_one, check_option_workflow, restricted_float
from ppanggolin.annotate import annotatePangenome, readAnnotations, getGeneSequencesFromFastas
from ppanggolin.cluster import clustering, readClustering
from ppanggolin.graph import computeNeighborsGraph
from ppanggolin.nem.rarefaction import makeRarefactionCurve
from ppanggolin.nem.partition import partition
from ppanggolin.formats import writePangenome, writeFlatFiles
from ppanggolin.figures import drawTilePlot, drawUCurve
from ppanggolin.info import printInfo


""" a global workflow that does everything in one go. """


def launch(args):
    check_option_workflow(args)
    pangenome = Pangenome()
    filename = mkFilename(args.basename, args.output, args.force)
    if args.anno:  # if the annotations are provided, we read from it
        readAnnotations(pangenome, args.anno, cpu=args.cpu, disable_bar=args.disable_prog_bar)
        writePangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is None:
            raise Exception("The gff/gbff provided did not have any sequence informations, "
                            "you did not provide clusters and you did not provide fasta file. "
                            "Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is not None:
            getGeneSequencesFromFastas(pangenome, args.fasta)

        if args.clusters is not None:
            readClustering(pangenome, args.clusters, disable_bar=args.disable_prog_bar)

        elif args.clusters is None:  # we should have the sequences here.
            clustering(pangenome, tmpdir=args.tmpdir, cpu=args.cpu, identity=args.identity, coverage=args.coverage,
                       mode=args.mode, defrag=not args.no_defrag, disable_bar=args.disable_prog_bar)
    elif args.fasta is not None:
        pangenome = Pangenome()
        annotatePangenome(pangenome, args.fasta, args.tmpdir, args.cpu, contig_filter=args.contig_filter,
                          disable_bar=args.disable_prog_bar)
        writePangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
        clustering(pangenome, tmpdir=args.tmpdir, cpu=args.cpu, identity=args.identity, coverage=args.coverage,
                   mode=args.mode, defrag=not args.no_defrag, disable_bar=args.disable_prog_bar)

    computeNeighborsGraph(pangenome, disable_bar=args.disable_prog_bar)

    partition(pangenome, tmpdir=args.tmpdir, cpu=args.cpu, kval=args.nb_of_partitions,
              disable_bar=args.disable_prog_bar)
    writePangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)

    if args.rarefaction:
        makeRarefactionCurve(pangenome, args.output, args.tmpdir, cpu=args.cpu, disable_bar=args.disable_prog_bar)
    if 1 < len(pangenome.organisms) < 5000:
        drawTilePlot(pangenome, args.output, nocloud=False if len(pangenome.organisms) < 500 else True)
    drawUCurve(pangenome, args.output)

    writeFlatFiles(pangenome, args.output, args.cpu, csv=True, genePA=True, gexf=True, light_gexf=True, projection=True,
                   json=True, stats=True, partitions=True)

    printInfo(filename, content=True)


def workflowSubparser(subparser):
    parser = subparser.add_parser("workflow", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Input arguments", description="The possible input arguments :")
    required.add_argument('--fasta', required=False, type=str,
                          help="A tab-separated file listing the organism names, "
                               "and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). "
                               "One line per organism. This option can be used alone.")
    required.add_argument('--anno', required=False, type=str,
                          help="A tab-separated file listing the organism names, and the gff filepath "
                               "of its annotations (the gffs can be compressed). One line per organism. "
                               "This option can be used alone IF the fasta sequences are in the gff files, "
                               "otherwise --fasta needs to be used.")
    required.add_argument("--clusters", required=False, type=str,
                          help="a tab-separated file listing the cluster names, the gene IDs, "
                               "and optionally whether they are a fragment or not.")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    optional.add_argument("--basename", required=False, default="pangenome", help="basename for the output file")
    optional.add_argument("--mode", required=False, default="1", choices=["0", "1", "2", "3"],
                          help="the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component),"
                               " 2: CD-HIT-like, 3: CD-HIT-like (lowmem)")
    optional.add_argument("--coverage", required=False, type=restricted_float, default=0.8,
                          help="Minimal coverage of the alignment for two proteins to be in the same cluster")
    optional.add_argument("--identity", required=False, type=restricted_float, default=0.8,
                          help="Minimal identity percent for two proteins to be in the same cluster")
    optional.add_argument("--rarefaction", required=False, action="store_true",
                          help="Use to compute the rarefaction curves (WARNING: can be time consuming)")
    optional.add_argument("-K", "--nb_of_partitions", required=False, default=-1, type=int,
                          help="Number of partitions to use. Must be at least 2. If under 2, "
                               "it will be detected automatically.")
    optional.add_argument("--defrag", required=False, action="store_true",
                          help=argparse.SUPPRESS)
    # This ensures compatibility with workflows built with the old option "defrag" when it was not the default
    optional.add_argument("--no_defrag", required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with their non-fragmented gene family.")
    optional.add_argument("--contig_filter", required=False, default=1, type=min_one, help=argparse.SUPPRESS)
    return parser
