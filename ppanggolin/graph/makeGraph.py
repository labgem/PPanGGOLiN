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
    logging.getLogger().debug(f"Ram used at the start : {getCurrentRAM()}")
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    computeNeighborsGraph(pangenome)
    writePangenome(pangenome, pangenome.file, args.force)


def graphSubparser(subparser):
    parser = subparser.add_parser("graph",help = "Create the pangenome graph")
    parser.add_argument("-p","--pangenome",type=str, help = "PPanGGOLiN binary pangenome files to read from and to add informations to. If provided the -o and --basename options will be ignored and data will be added to this pangenome instead.")
    # parser.add_argument('-org', '--organisms', type=str, help="""
    # File: A tab-delimited file containing at least 2 mandatory fields per row and as many optional fields as the number of circular contigs. 
    # Each row corresponds to an organism to be added to the pangenome.
    # The first field is the organism ID.
    # The organism ID can be any string but must be unique and can't contain any space, quote, double quote
    # The second field is the gff file containing the annotations associated to the organism. 
    # This path can be absolute or relative. 
    # The gff file must contain an ID for each line.
    # Accepted types are CDS or xRNA (rRNA,tRNA,tmRNA) in the type column.
    # The contig ID and gene ID can be any string but must be unique and can't contain any space, quote, double quote, pipe
    # (optional): The next fields contain the name of perfectly assembled circular contigs (to take in account the link between the first and the last gene in the graph). 
    # In this case, it is mandatory to provide the contig size in the gff files either by adding a "region" type having the correct contig ID attribute or using a '##sequence-region' pragma.
    # """)
    # parser.add_argument('-gf', '--gene_families', type=str, help="""
    # File: A tab-delimited file containing the gene families. Each row contains 2 or 3 fields.
    # The first field is the family ID. 
    # The second field is the gene IDs associated with this family. 
    # The third field (optional) is a flag "F" to specify if the gene is a gene fragment (empty otherwise).
    # If several consecutive genes belonging to the same gene families have the flag, then the no reflexive links are drawn between them in the graph.
    # The family ID can be any string but must be unique and can't contain any space, quote, double quote and reserved word. 
    # The family ID is often the name of the most representative sequence of the family.
    # Gene IDs can be any string corresponding to the IDs of the gff files. They must be uniques and can't contain any spaces, quote, double quote, pipe and reserved words.
    # """)
    parser.add_argument('-r', '--remove_high_copy_number_families', type=int, default=0, help="""
    Positive Number: Remove families having a number of copy of gene in a single organism above or equal to this threshold in at least one organism (0 or negative values are ignored). 
    """)
    parser.add_argument('-s', '--infer_singletons', default=False, action="store_true", help="""
    Flag: If a gene id found in a gff file is absent of the gene families file, a gene families with a single gene will be created.
    if this argument is not set, the program will raise a KeyError exception if a gene id found in a gff file is absent of the gene families file.""")
    return parser