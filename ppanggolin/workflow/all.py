#!/usr/bin/env python3
#coding:utf-8

#default libraries
import os
import time
import argparse
import logging

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mkFilename, min_one
from ppanggolin.annotate import annotatePangenome, readAnnotations, getGeneSequencesFromFastas
from ppanggolin.cluster import clustering, readClustering
from ppanggolin.graph import computeNeighborsGraph
from ppanggolin.nem.rarefaction import makeRarefactionCurve
from ppanggolin.nem.partition import partition
from ppanggolin.formats import writePangenome, writeFlatFiles
from ppanggolin.figures import drawTilePlot, drawUCurve
from ppanggolin.info import printInfo
from ppanggolin.RGP.genomicIsland import predictRGP
from ppanggolin.RGP.spot import predictHotspots
from ppanggolin.mod import predictModules
### a global workflow that does everything in one go.

def launch(args):
    pangenome = Pangenome()
    filename = mkFilename(args.basename, args.output, args.force)
    if args.anno:#if the annotations are provided, we read from it
        getSeq = True
        if args.clusters is not None:
            getSeq = False
        start_anno = time.time()
        readAnnotations(pangenome, args.anno, cpu = args.cpu, getSeq = getSeq, show_bar=args.show_prog_bars)
        annotime = time.time() - start_anno
        start_writing = time.time()
        writePangenome(pangenome, filename, args.force, show_bar=args.show_prog_bars)
        writing_time = time.time() - start_writing
        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is None:
            raise Exception("The gff/gbff provided did not have any sequence informations, you did not provide clusters and you did not provide fasta file. Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is not None:
            getGeneSequencesFromFastas(pangenome, args.fasta)
        start_clust = time.time()
        if args.clusters is not None:
            readClustering(pangenome, args.clusters, show_bar=args.show_prog_bars)

        elif args.clusters is None:#we should have the sequences here.
            clustering(pangenome, args.tmpdir, args.cpu, defrag=not args.no_defrag, show_bar=args.show_prog_bars)
        clust_time = time.time() - start_clust
    elif args.fasta is not None:
        start_anno = time.time()
        annotatePangenome(pangenome, args.fasta, args.tmpdir, args.cpu, contig_filter=args.contig_filter, show_bar=args.show_prog_bars)
        annotime = time.time() - start_anno
        start_writing = time.time()
        writePangenome(pangenome, filename, args.force, show_bar=args.show_prog_bars)
        writing_time = time.time() - start_writing
        start_clust = time.time()
        clustering(pangenome, args.tmpdir, args.cpu, defrag=not args.no_defrag, show_bar=args.show_prog_bars)
        clust_time = time.time() - start_clust

    writePangenome(pangenome, filename, args.force, show_bar=args.show_prog_bars)
    start_graph = time.time()
    computeNeighborsGraph(pangenome, show_bar=args.show_prog_bars)
    graph_time = time.time() - start_graph

    start_part = time.time()
    partition(pangenome, tmpdir = args.tmpdir, cpu = args.cpu, K=args.nb_of_partitions, show_bar=args.show_prog_bars)
    part_time = time.time() - start_part

    start_writing = time.time()
    writePangenome(pangenome, filename, args.force, show_bar=args.show_prog_bars)
    writing_time = writing_time + time.time() - start_writing

    start_regions = time.time()
    predictRGP(pangenome, show_bar=args.show_prog_bars)
    regions_time = time.time() - start_regions

    start_spots = time.time()
    predictHotspots(pangenome, args.output, show_bar=args.show_prog_bars)
    spot_time = time.time() - start_spots

    start_mods = time.time()
    predictModules(pangenome = pangenome, cpu = args.cpu, tmpdir = args.tmpdir, show_bar=args.show_prog_bars)
    mod_time = time.time() - start_mods

    start_writing = time.time()
    writePangenome(pangenome, filename, args.force, show_bar=args.show_prog_bars)
    writing_time = writing_time + time.time() - start_writing

    if not args.only_pangenome:
        if args.rarefaction:
            makeRarefactionCurve(pangenome,args.output, args.tmpdir, cpu=args.cpu, show_bar=args.show_prog_bars)
        if len(pangenome.organisms) > 1 and len(pangenome.organisms) < 5000:
            drawTilePlot(pangenome, args.output, nocloud = False if len(pangenome.organisms) < 500 else True)
        drawUCurve(pangenome, args.output)

        start_desc = time.time()
        writeFlatFiles(pangenome, args.output, args.cpu, csv = True, genePA=True, gexf=True, light_gexf = True, projection=True, json = True, stats = True, partitions = True, regions = True, spots=True, borders=True, spot_modules=True, modules=True)
        desc_time = time.time() - start_desc

    logging.getLogger().info(f"Annotation took : {round(annotime,2)} seconds")
    logging.getLogger().info(f"Clustering took : {round(clust_time,2)} seconds")
    logging.getLogger().info(f"Building the graph took : {round(graph_time,2)} seconds")
    logging.getLogger().info(f"Partitionning the pangenome took : {round(part_time,2)} seconds")
    logging.getLogger().info(f"Predicting RGP took : {round(regions_time,2)} seconds")
    logging.getLogger().info(f"Gathering RGP into spots took : {round(spot_time,2)} seconds")
    logging.getLogger().info(f"Predicting modules took : {round(start_mods, 2)} seconds")
    logging.getLogger().info(f"Writing the pangenome data in HDF5 took : {round(writing_time,2)} seconds")
    printInfo(filename, content = True)


def allSubparser(subparser):
    parser = subparser.add_parser("all", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Input arguments", description = "The possible input arguments :")
    required.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism. This option can be used alone.")
    required.add_argument('--anno', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used.")
    required.add_argument("--clusters",required=False, type=str, help = "a tab-separated file listing the cluster names, the gene IDs, and optionnally whether they are a fragment or not.")

    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--basename",required = False, default = "pangenome", help = "basename for the output file")
    optional.add_argument("-K","--nb_of_partitions",required=False, default=-1, type=int, help = "Number of partitions to use. Must be at least 3. If under 3, it will be detected automatically.")
    optional.add_argument("--defrag", required=False, action = "store_true", help = argparse.SUPPRESS)##This ensures compatibility with workflows built with the old option "defrag" when it was not the default
    optional.add_argument("--no_defrag",required=False, action="store_true", help = "DO NOT Realign gene families to link fragments with their non-fragmented gene family.")
    optional.add_argument("--only_pangenome", required=False, action="store_true", help = "Only generate the HDF5 pangenome file")
    optional.add_argument("--contig_filter",required=False, default=1, type=min_one, help = "remove contigs that are smaller than this length")
    return parser