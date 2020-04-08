#!/usr/bin/env python3
#coding:utf-8

#default libraries
import os
import time
import argparse
import logging

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
from ppanggolin.RGP.genomicIsland import predictRGP
from ppanggolin.RGP.spot import predictHotspots
### a global workflow that does everything in one go.

def launch(args):
    pangenome = Pangenome()
    filename = mkFilename(args.basename, args.output, args.force)
    if args.anno:#if the annotations are provided, we read from it
        getSeq = True
        if args.clusters is not None:
            getSeq = False
        start_anno = time.time()
        readAnnotations(pangenome, args.anno, cpu = args.cpu, getSeq = getSeq)
        annotime = time.time() - start_anno
        start_writing = time.time()
        writePangenome(pangenome, filename, args.force)
        writing_time = time.time() - start_writing
        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is None:
            raise Exception("The gff/gbff provided did not have any sequence informations, you did not provide clusters and you did not provide fasta file. Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is not None:
            getGeneSequencesFromFastas(pangenome, args.fasta)
        start_clust = time.time()
        if args.clusters is not None:
            readClustering(pangenome, args.clusters)

        elif args.clusters is None:#we should have the sequences here.
            clustering(pangenome, args.tmpdir, args.cpu, defrag=args.defrag)
        clust_time = time.time() - start_clust
    elif args.fasta is not None:
        start_anno = time.time()
        annotatePangenome(pangenome, args.fasta, args.tmpdir, args.cpu)
        annotime = time.time() - start_anno
        start_writing = time.time()
        writePangenome(pangenome, filename, args.force)
        writing_time = time.time() - start_writing
        start_clust = time.time()
        clustering(pangenome, args.tmpdir, args.cpu, defrag=args.defrag)
        clust_time = time.time() - start_clust

    writePangenome(pangenome, filename, args.force)
    start_graph = time.time()
    computeNeighborsGraph(pangenome)
    graph_time = time.time() - start_graph

    start_part = time.time()
    partition(pangenome, tmpdir = args.tmpdir, cpu = args.cpu, K=args.nb_of_partitions)
    part_time = time.time() - start_part

    start_writing = time.time()
    writePangenome(pangenome, filename, args.force)
    writing_time = writing_time + time.time() - start_writing

    start_regions = time.time()
    predictRGP(pangenome)
    regions_time = time.time() - start_regions

    start_spots = time.time()
    predictHotspots(pangenome, args.output, interest=args.interest)
    spot_time = time.time() - start_spots

    start_writing = time.time()
    writePangenome(pangenome, filename, args.force)
    writing_time = writing_time + time.time() - start_writing

    if args.rarefaction:
        makeRarefactionCurve(pangenome,args.output, args.tmpdir, cpu=args.cpu)
    if len(pangenome.organisms) > 1 and len(pangenome.organisms) < 5000:
        drawTilePlot(pangenome, args.output, nocloud = False if len(pangenome.organisms) < 500 else True)
    drawUCurve(pangenome, args.output)

    start_desc = time.time()
    writeFlatFiles(pangenome, args.output, args.cpu, csv = True, genePA=True, gexf=True, light_gexf = True, projection=True, json = True, stats = True, partitions = True, regions = True, spots=True)
    desc_time = time.time() - start_desc

    logging.getLogger().info(f"Annotation took : {round(annotime,2)} seconds")
    logging.getLogger().info(f"Clustering took : {round(clust_time,2)} seconds")
    logging.getLogger().info(f"Building the graph took : {round(graph_time,2)} seconds")
    logging.getLogger().info(f"Partitionning the pangenome took : {round(part_time,2)} seconds")
    logging.getLogger().info(f"Predicting RGP took : {round(regions_time,2)} seconds")
    logging.getLogger().info(f"Gathering RGP into spots took : {round(spot_time,2)} seconds")
    logging.getLogger().info(f"Writing the pangenome data in HDF5 took : {round(writing_time,2)} seconds")
    logging.getLogger().info(f"Writing descriptive files for the pangenome took : {round(desc_time,2)} seconds")
    printInfo(filename, content = True)


def panRGPSubparser(subparser):
    parser = subparser.add_parser("panrgp", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Input arguments", description = "The possible input arguments :")
    required.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism. This option can be used alone.")
    required.add_argument('--anno', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used.")
    required.add_argument("--clusters",required=False, type=str, help = "a tab-separated file listing the cluster names, the gene IDs, and optionnally whether they are a fragment or not.")

    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--basename",required = False, default = "pangenome", help = "basename for the output file")
    optional.add_argument("--rarefaction", required=False, action = "store_true", help = "Use to compute the rarefaction curves (WARNING: can be time consumming)")
    optional.add_argument("-K","--nb_of_partitions",required=False, default=-1, type=int, help = "Number of partitions to use. Must be at least 3. If under 3, it will be detected automatically.")
    optional.add_argument("--interest",required=False, type=str, default="",help = "Comma separated list of elements to flag when drawing and writing hotspots")
    optional.add_argument("--defrag",required=False, action="store_true", help = "Realign gene families to associated fragments with their non-fragmented gene family.")

    return parser