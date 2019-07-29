#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import logging
import os

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import getCurrentRAM
from .makeNeighborsGraph import readClustering, readAnnotations, computeNeighborsGraph, linkGeneFamsToGenes


def launchFromDB(args):
    """
        :param datapath: a path to PPanGGOLiN datafiles containing a previously computed pangenome, either with annotations only, or annotations and associated gene families.
        :type file:
    """
    pangenome = Pangenome()
    logging.getLogger().debug(f"RAM usage prior starting : {getCurrentRAM()}")
    pangFiles = checkDatapath(args.pangenome)#get which db are provided
    #if graph is there already, exit as this step is useless.
    if "graph" in pangFiles:
        logging.getLogger().error(f"The graph was already computed and is in the provided pangenome file(s) : '{args.pangenome}'. Exiting.")
        exit(1)
    if len(pangFiles) == 2:
        raise FileNotFoundError(f"No files were found under the name {args.pangenome}.")
    logging.getLogger().info("Reading ppanggolin binary file(s)")
    #we want to read annotations first as what will be read will impact the gene families reading
    if "annotations" in pangFiles and "families" in pangFiles:
        readAnnotationsBinary(pangFiles["annotations"], pangenome)
        logging.getLogger().debug(f"RAM after loading annotations: {getCurrentRAM()}")
        if len(pangenome.geneFamilies) == 0:#there was no gene families in the gene annotations
            linkGeneFamsToGenes(readGeneFamiliesBinaryClustering(pangFiles["families"], pangenome), pangenome)
        else:
            readGeneFamiliesBinary(pangFiles["families"], pangenome)
        logging.getLogger().info(f"Linking gene families together into a neighbors graph...")
        computeNeighborsGraph(pangenome)
        logging.getLogger().debug(f"RAM with neighbors graph : {getCurrentRAM()}")
        logging.getLogger().info("Done creating the graph")
        #write only the graph.
        writeGraphBinary(pangenome,pangFiles["directory"], pangFiles["basename"])
    elif "annotations" in pangFiles and args.gene_families is not None:
        readAnnotationsBinary(pangFiles["annotations"], pangenome)
        logging.getLogger().debug(f"RAM after loading annotations: {getCurrentRAM()}")
        linkGeneFamsToGenes(readClustering(pangenome, args.gene_families), pangenome)
        logging.getLogger().info(f"Linking gene families together into a neighbors graph...")
        computeNeighborsGraph(pangenome)
        logging.getLogger().debug(f"RAM with neighbors graph : {getCurrentRAM()}")
        logging.getLogger().info("Done creating the graph")
        #rewrite everything as there are new informations also in 'annotations'
        writeBinaries(pangenome, pangFiles["directory"], pangFiles["basename"])
    elif "families" in pangFiles and args.organisms is not None:
        logging.getLogger().info(f"Loading annotations and gene families...")
        readAnnotations(pangenome, args.organisms, args.remove_high_copy_number_families, args.infer_singletons, args.add_rna_to_the_pangenome, readGeneFamiliesBinaryClustering(pangFiles["families"], pangenome))
        logging.getLogger().info("Done reading annotation files")
        logging.getLogger().info(f"Linking gene families together into a neighbors graph...")
        computeNeighborsGraph(pangenome)
        logging.getLogger().debug(f"RAM with neighbors graph : {getCurrentRAM()}")
        logging.getLogger().info("Done creating the graph")
        #write everything as things in the gene families may have changed
        writeBinaries(pangenome, pangFiles["directory"], pangFiles["basename"])
    else:#combination of option is impossible.
        pass

def launchFromFiles(organisms_file, families_tsv_file, lim_occurence, infer_singletons, add_rna_to_the_pangenome, output, basename):
    """
        :param organisms_file: a file listing organims by compute, first column is organism name, second is path to gff file and optionnally other other to provide the name of circular contig
        :param families_tsv_file: a file listing families. The first element is the family identifier (by convention, we advice to use the identifier of the average gene of the family) and then the next elements are the identifiers of the genes belonging to this family.
        :param lim_occurence: a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the families_repeted attribute.
        :param add_rna_to_the_pangenome: a bool specifying if the rna genes must be added to the pangenome or not.
        :param infer_singletons: a bool specifying if singleton must be explicitely present in the families_tsv_file (False) or if single gene in gff files must be automatically infered as a singleton family (True)
        :param output: a directory where to store all the pangenome files
        :param basename: a basename for the pangenome files
        :type file: 
        :type file: 
        :type int: 
        :type bool: 
        :type int:
        :type bool: 
        :type directory:
        :type str:
    """
    pangenome = Pangenome()
    logging.getLogger().debug(f"RAM usage prior starting : {getCurrentRAM()}")
    logging.getLogger().info(f"Loading annotations and gene families...")
    readAnnotations(pangenome, organisms_file, lim_occurence, infer_singletons, add_rna_to_the_pangenome, readClustering(pangenome, families_tsv_file))
    logging.getLogger().info("Done reading annotation files")
    logging.getLogger().info(f"Linking gene families together into a neighbors graph...")
    computeNeighborsGraph(pangenome)
    logging.getLogger().debug(f"RAM with neighbors graph : {getCurrentRAM()}")
    logging.getLogger().info("Done creating the graph")
    writeBinaries(pangenome, output, basename)

def launch(args):
    """
        This function is called only when using the module 'create' from the command line.
    """
    #check if we were provided with an existing directory or not, and reacting accordingly
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    elif not args.force:
        logging.getLogger().error(args.output+" already exists")
        exit(1)
    #check input type (file or db)
    #DB (annotation DB and cluster DB ?)
    if args.pangenome is not None:
        launchFromDB(args)
    #Files only
    elif args.organisms is not None and args.gene_families is not None:
        launchFromFiles(args.organisms, args.gene_families, args.remove_high_copy_number_families, args.infer_singletons, args.add_rna_to_the_pangenome, args.output, args.basename)
        # write all binaries (annotations, graph, gene families)
        

def createSubparser(subparser):
    parser = subparser.add_parser("create",help = "Create the pangenome graph")
    parser.add_argument("-p","--pangenome",type=str, help = "PPanGGOLiN binary pangenome files to read from and to add informations to. If provided the -o and --basename options will be ignored and data will be added to this pangenome instead.")
    parser.add_argument('-org', '--organisms', type=str, help="""
    File: A tab-delimited file containing at least 2 mandatory fields per row and as many optional fields as the number of circular contigs. 
    Each row corresponds to an organism to be added to the pangenome.
    The first field is the organism ID.
    The organism ID can be any string but must be unique and can't contain any space, quote, double quote
    The second field is the gff file containing the annotations associated to the organism. 
    This path can be absolute or relative. 
    The gff file must contain an ID for each line.
    Accepted types are CDS or xRNA (rRNA,tRNA,tmRNA) in the type column.
    The contig ID and gene ID can be any string but must be unique and can't contain any space, quote, double quote, pipe
    (optional): The next fields contain the name of perfectly assembled circular contigs (to take in account the link between the first and the last gene in the graph). 
    In this case, it is mandatory to provide the contig size in the gff files either by adding a "region" type having the correct contig ID attribute or using a '##sequence-region' pragma.
    """)
    parser.add_argument('-gf', '--gene_families', type=str, help="""
    File: A tab-delimited file containing the gene families. Each row contains 2 or 3 fields.
    The first field is the family ID. 
    The second field is the gene IDs associated with this family. 
    The third field (optional) is a flag "F" to specify if the gene is a gene fragment (empty otherwise).
    If several consecutive genes belonging to the same gene families have the flag, then the no reflexive links are drawn between them in the graph.
    The family ID can be any string but must be unique and can't contain any space, quote, double quote and reserved word. 
    The family ID is often the name of the most representative sequence of the family.
    Gene IDs can be any string corresponding to the IDs of the gff files. They must be uniques and can't contain any spaces, quote, double quote, pipe and reserved words.
    """)
    parser.add_argument('-r', '--remove_high_copy_number_families', type=int, default=0, help="""
    Positive Number: Remove families having a number of copy of gene in a single family above or equal to this threshold in at least one organism (0 or negative values are ignored). 
    """)
    parser.add_argument('-s', '--infer_singletons', default=False, action="store_true", help="""
    Flag: If a gene id found in a gff file is absent of the gene families file, a gene families with a single gene will be created.
    if this argument is not set, the program will raise a KeyError exception if a gene id found in a gff file is absent of the gene families file.""")
    parser.add_argument("-ra", "--add_rna_to_the_pangenome", default = False, action="store_true", help = """
    Flag: If the specified the xRNA (rRNA,tRNA,tmRNA...) are added to the pangenome graph, (like proteins, RNA genes need to be clustered)""")
    parser.add_argument("--basename", required=False, type=str, default = "pangenomeGraph", help = "Used to name output files")
    return parser