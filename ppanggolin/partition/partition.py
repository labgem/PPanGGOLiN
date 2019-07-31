#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import logging


#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import getCurrentRAM
from ppanggolin.formats import readPangenome, writePangenome
#cython library
import nem_stats

def evaluate_nb_partitions(pangenome, beta, sm_degree, free_dispersion, chunk_size, Q, Qrange, ICL_margin, cpu, tmpdir):
    ChosenQne = 3
    logging.getLogger().info("Estimating the optimal number of partitions...")





    return ChosenQne

def checkPangenomePartition(pangenome):
    if pangenome.status["genomesAnnotated"] in ["Computed","Loaded"]:
        pass
    elif pangenome.status["genomesAnnotated"] == "inFile":
        readPangenome(pangenome, annotation = True)
    else:
        raise Exception("You want to partition an unannotated pangenome")
    if pangenome.status["genesClustered"] in ["Computed","Loaded"]:
        pass
    elif pangenome.status["genesClustered"] == "inFile":
        readPangenome(pangenome, geneFamilies= True)
    else:
        raise Exception("You want to partition a pangenome whose genes have not been clustered")
    if pangenome.status["neighborsGraph"] in ["Computed","Loaded"]:
        pass
    elif pangenome.status["neighborsGraph"] == "inFile":
        readPangenome(pangenome, graph=True)#whether it is faster to compute it or to load it will have to be checked on bigger graphs.
        #also maybe it's not even that useful since we might recompute it for all subpangenomes?
    else:
        raise Exception("You want to partition a pangenome whose neighbors graph has not been computed.")

def partition(pangenome, beta, sm_degree, free_dispersion, chunk_size, Q, Qrange, ICL_margin, cpu, tmpdir):
    checkPangenomePartition(pangenome)

    if Q < 3:
        Q = evaluate_nb_partitions(pangenome, beta, sm_degree, free_dispersion, chunk_size, Q, Qrange, ICL_margin, cpu, tmpdir)
    
    
    pass

def launch(args):
    """
        main code when launch partition from the command line.
    """
    logging.getLogger().debug(f"Ram used at the start : {getCurrentRAM()}")
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    partition(pangenome, args.beta, args.max_degree_smoothing, args.free_dispersion, args.chunk_size, args.nb_of_partitions, args.qrange, args.ICL_margin, args.cpu, args.tmpdir)



def partitionSubparser(subparser):
    parser = subparser.add_parser("partition",help = "Partition the pangenome graph")
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument("-b","--beta", required = False, default = 2.5, type = float, help = "beta is the strength of the smoothing using the graph topology during partitionning. 0 will deactivate spatial smoothing.")
    optional.add_argument("-ms","--max_degree_smoothing",required = False, default = float("inf"), help = "max. degree of the nodes to be included in the smoothing process.")
    #soft core ???
    optional.add_argument("-fd","--free_dispersion",required = False, default = False, action = "store_true",help = "use if the dispersion around the centroid vector of each partition during must be free. It will be the same for all organisms by default.")
    optional.add_argument("-ck","--chunk_size",required=False, default = 500, type = int, help = "Size of the chunks when performing partitionning using chunks of organisms. Chunk partitionning will be used automatically if the number of genomes is above this number.")
    optional.add_argument("-Q","--nb_of_partitions",required=False, default=-1, type=int, help = "Number of partitions to use. Must be at least 3. If under 3, it will be detected automatically.")
    optional.add_argument("-Qmm","--qrange",nargs=2,required = False, type=int, default=[3,20], help="Range of Q values to test when detecting Q automatically. Default between 3 and 20.")
    optional.add_argument("-im","--ICL_margin",required = False, type = float, default = 0.05, help = "Q is detected automatically by maximizing ICL. However at some point the ICL reaches a plateau. Therefore we are looking for the minimal value of Q without significative gain from the larger values of Q measured by ICL. For that we take the lowest Q that is found within a given 'margin' of the maximal ICL value. Basically, change this option only if you truly understand it, otherwise just leave it be.")
    optional.add_argument("--draw_ICL", required =False, default = False, action="store_true",help = "Use if you can to draw the ICL curve for all of the tested Q values. Will not be done if Q is given.")
    #use former partitionning ?????
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism.")
    
    return parser