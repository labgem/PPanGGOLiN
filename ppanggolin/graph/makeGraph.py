#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import logging
import os
#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import getCurrentRAM
from ppanggolin.formats import readPangenome, writePangenome

def checkPangenomeForNeighborsGraph(pangenome):
    """
        Checks the pangenome for neighbors graph computing.
    """
    if pangenome.status["genomesAnnotated"] in ["Computed","Loaded"] and pangenome.status["genesClustered"] in ["Computed","Loaded"]:
        pass#nothing to do, can just continue.
    elif pangenome.status["genomesAnnotated"] == "inFile" and pangenome.status["genesClustered"] == "inFile":
        readPangenome(pangenome, annotation = True, geneFamilies=True)
    else:
        #You probably can use readPangenome anyways.
        msg = "Dev : You are probably writing a new workflow with a combination that I did not test. You can probably use readPangenome instead of raising this Error. However please test it carefully.\n"
        msg+=" User : I have no idea how you got there. You probably did something unexpected. Contact the devs (or try something else).\n"
        raise NotImplementedError(msg)

def computeNeighborsGraph(pangenome):
    """
        Creates the Pangenome Graph. Will either load the informations from the pangenome file if they are not loaded, or use the informations loaded if they are.
    """
    checkPangenomeForNeighborsGraph(pangenome)
    logging.getLogger().info("Computing the neighbors graph...")
    bar = tqdm(pangenome.organisms, total = len(pangenome.organisms), unit = "organism")
    for org in bar:
        bar.set_description(f"Processing {org.name}")
        bar.refresh()
        for contig in org.contigs:
            prev = None
            for gene in contig.genes:
                if prev is not None and not gene.family.removed:
                    if not (prev.family == gene.family and (prev.is_fragment or gene.is_fragment)):
                        pangenome.addEdge(gene, prev)
                prev = gene
            if contig.is_circular:
                pangenome.addEdge(contig.genes[0],prev)
    logging.getLogger().info("Done making the neighbors graph.")
    pangenome.status["neighborsGraph"] = "Computed"


def launch(args):
    if args.remove_high_copy_number_families != 0:
        raise NotImplementedError()
    logging.getLogger().debug(f"Ram used at the start : {getCurrentRAM()}")
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    computeNeighborsGraph(pangenome)
    writePangenome(pangenome, pangenome.file, args.force)


def graphSubparser(subparser):
    parser = subparser.add_parser("graph",help = "Create the pangenome graph")
    parser.add_argument("-p","--pangenome",type=str, help = "PPanGGOLiN binary pangenome files to read from and to add informations to. If provided the -o and --basename options will be ignored and data will be added to this pangenome instead.")
    parser.add_argument('-r', '--remove_high_copy_number_families', type=int, default=0, help="""
    Positive Number: Remove families having a number of copy of gene in a single organism above or equal to this threshold in at least one organism (0 or negative values are ignored). 
    """)
    return parser