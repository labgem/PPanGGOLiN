#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging

#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import getCurrentRAM
from ppanggolin.formats import readPangenome, writePangenome, ErasePangenome

def checkPangenomeFormerGraph(pangenome, force):
    """ checks pangenome status and .h5 files for former neighbors graph, delete it if allowed or raise an error """
    if pangenome.status["neighborsGraph"] == "inFile" and force == False:
        raise Exception("You are trying to make a neighbors graph that is already built. If you REALLY want to do that, use --force (it will erase everything except annotation data !)")
    elif pangenome.status["neighborsGraph"] == "inFile" and force == True:
        ErasePangenome(pangenome, graph = True)

def checkPangenomeForNeighborsGraph(pangenome, force):
    """
        Checks the pangenome for neighbors graph computing.
    """
    checkPangenomeFormerGraph(pangenome, force)
    if pangenome.status["genomesAnnotated"] in ["Computed","Loaded"] and pangenome.status["genesClustered"] in ["Computed","Loaded"]:
        pass#nothing to do, can just continue.
    elif pangenome.status["genomesAnnotated"] == "inFile" and pangenome.status["genesClustered"] == "inFile":
        readPangenome(pangenome, annotation = True, geneFamilies=True)
    else:
        #You probably can use readPangenome anyways.
        msg = "Dev : You are probably writing a new workflow with a combination that I did not test. You can probably use readPangenome instead of raising this Error. However please test it carefully.\n"
        msg+=" User : I have no idea how you got there. You probably did something unexpected. Contact the devs (or try something else).\n"
        raise NotImplementedError(msg)

def remove_high_copy_number(pangenome, number):
    """ removes families present more than 'number' times from the pangenome graph"""
    for fam in pangenome.geneFamilies:
        for gene_list in fam.getOrgDict().values():
            if len(gene_list) >= number:
                fam.removed = True


def computeNeighborsGraph(pangenome, remove_copy_number = 0, force = False):
    """
        Creates the Pangenome Graph. Will either load the informations from the pangenome file if they are not loaded, or use the informations loaded if they are.
    """
    checkPangenomeForNeighborsGraph(pangenome, force)

    if remove_copy_number > 0 :
        remove_high_copy_number(pangenome, remove_copy_number)

    logging.getLogger().info("Computing the neighbors graph...")
    bar = tqdm(pangenome.organisms, total = len(pangenome.organisms), unit = "organism")
    for org in bar:
        bar.set_description(f"Processing {org.name}")
        bar.refresh()
        for contig in org.contigs:
            prev = None
            for gene in contig.genes:
                if not gene.family.removed:
                    if prev is not None:
                        if not (prev.family == gene.family and (prev.is_fragment or gene.is_fragment)):
                            pangenome.addEdge(gene, prev)
                    prev = gene
            if contig.is_circular and len(contig.genes) > 0:
                pangenome.addEdge(contig.genes[0],prev)
    logging.getLogger().info("Done making the neighbors graph.")
    pangenome.status["neighborsGraph"] = "Computed"

    pangenome.parameters["graph"] = {}
    pangenome.parameters["graph"]["removed_high_copy_number_families"] = False
    if remove_copy_number > 0:
        pangenome.parameters["graph"]["removed_high_copy_number_families"] = True
        pangenome.parameters["graph"]["removed_high_copy_number_of_families_above"] = remove_copy_number


def launch(args):
    logging.getLogger().debug(f"Ram used at the start : {getCurrentRAM()}")
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    computeNeighborsGraph(pangenome, args.remove_high_copy_number, args.force)
    writePangenome(pangenome, pangenome.file, args.force)


def graphSubparser(subparser):
    parser = subparser.add_parser("graph")
    parser.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    parser.add_argument('-r', '--remove_high_copy_number', type=int, default=0, help="""Positive Number: Remove families having a number of copy of gene in a single organism above or equal to this threshold in at least one organism (0 or negative values are ignored).""")
    return parser