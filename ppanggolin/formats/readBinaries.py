#!/usr/bin/env python3
#coding:utf-8

#default libraries
from pathlib import Path
import logging
import sys 
from multiprocessing import Pool

#installed libraries
from tqdm import tqdm
import tables

#local libraries
from ppanggolin.utils import getCurrentRAM, read_compressed_or_not
from ppanggolin.pangenome import Pangenome, GeneFamily, Edge
from ppanggolin.genome import Organism, Contig, Gene


def getStatus(pangenome, pangenomeFile):
    """
        Checks which elements are already present in the file.
    """
    h5f = tables.open_file(pangenomeFile,"r")
    logging.getLogger().info("Getting the current pangenome's status")
    statusGroup = h5f.root.status
    if statusGroup._v_attrs.genomesAnnotated:
        pangenome.status["genomesAnnotated"] = "inFile"
    if statusGroup._v_attrs.genesClustered:
        pangenome.status["genesClustered"] = "inFile"
    if statusGroup._v_attrs.geneSequences:
        pangenome.status["geneSequences"] = "inFile"
    if statusGroup._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "inFile"
    if statusGroup._v_attrs.NeighborsGraph:
        pangenome.status["NeighborsGraph"] = "inFile"
    if statusGroup._v_attrs.Partitionned:
        pangenome.status["Partitionned"] = "inFile"
    h5f.close()

def read_chunks(table, chunk=10000):
    """
        Reading entirely the provided table chunk per chunk to limit RAM usage.
    """
    for i in range(0,  table.nrows, chunk):
        for row in table.read(start = i, stop = i + chunk):
            yield row

def getGeneSequencesFromFile(pangenome, fileObj):
    """
        Writes the CDS sequences of the Pangenome object to a tmpFile object
        Loads the sequences from a .h5 pangenome file
    """
    logging.getLogger().info("Extracting and writing all of the CDS sequences from a .h5 pangenome file to a fasta file")
    h5f = tables.open_file(pangenome.file,"r", driver_core_backing_store=0)
    table = h5f.root.geneSequences
    bar =  tqdm(range(table.nrows), unit="gene")
    for row in read_chunks(table):#reading the table chunk per chunk otherwise RAM dies on big pangenomes
        if row[2] == b"CDS":
            fileObj.write('>' + row[1].decode() + "\n")
            fileObj.write(row[0].decode() + "\n")
        bar.update()
    fileObj.flush()
    bar.close()
    h5f.close()

def launchReadOrganism(args):
    return readOrganism(*args)

def readOrganism(pangenome, orgName, contigDict, link = False):
    org = Organism(orgName)
    for contigName, geneList in contigDict.items():
        contig = org.addContig(contigName, is_circular=geneList[0][0][0])
        for row in geneList:
            if link:#if the gene families are already computed/loaded the gene exists.
                gene = pangenome.getGene(row[1][0].decode())
            else:#else creating the gene.
                gene = Gene(row[1][0].decode())
            gene.fill_annotations(
                start = row[1][6],
                stop =row[1][7],
                strand =  row[1][8].decode(),
                geneType = row[1][9].decode(),
                position = row[1][4],
                genetic_code=row[1][1],
                name = row[1][3].decode(),
                product = row[1][5].decode())
            gene.is_fragment = row[1][2]
            gene.fill_parents(org, contig)
            if gene.type == "CDS":
                contig.addGene(gene)
            elif "RNA" in gene.type.upper():
                contig.addRNA(gene)
            else:
                raise Exception(f"A strange type ({gene.type}), which we do not know what to do with, was met.")
    pangenome.addOrganism(org)

def readGraph(pangenome, h5f):
    raise NotImplementedError()

def readGeneFamilies(pangenome, h5f):
    table = h5f.root.geneFamilies

    link = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] else False

    bar = tqdm(range(table.nrows), unit = "gene")
    for row in read_chunks(table):
        fam = pangenome.addGeneFamily(row[1].decode())
        if link:#linking if we have loaded the annotations
            geneObj = pangenome.getGene(row[0].decode())
        else:#else, no
            geneObj = Gene(row[1].decode())
        fam.addGene(geneObj)
        bar.update()
    bar.close()
    pangenome.status["genesClustered"] = "Loaded"

def readAnnotation(pangenome, h5f, filename):
    annotations = h5f.root.annotations
    
    table = annotations.genes
    bar = tqdm(range(table.nrows), unit="gene")
    pangenomeDict = {}
    for row in read_chunks(table):
        try:
            pangenomeDict[row[2].decode()][row[0][1].decode()].append(row)#new gene, seen contig, seen org
        except:
            try:
                pangenomeDict[row[2].decode()][row[0][1].decode()] = [row]#new contig, seen org
            except:
                pangenomeDict[sys.intern(row[2].decode())] = { row[0][1].decode() : [row]}#new org
        bar.update()
    bar.close()

    link = True if pangenome.status["genesClustered"] in ["Computed","Loaded"] else False

    bar = tqdm(range(len(pangenomeDict)), unit = "organism")
    for orgName, contigDict in pangenomeDict.items():
        readOrganism(pangenome, orgName, contigDict, link)
        bar.update()
    bar.close()
    pangenome.status["genomesAnnotated"] = "Loaded"

def readPangenome(pangenome, annotation = False, geneFamilies = False, graph = False):
    """
        Reads a previously written pangenome, with all of its parts.
    """
    # compressionFilter = tables.Filters(complevel=1, complib='blosc:lz4')
    if hasattr(pangenome,"file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename,"r")
    if annotation:
        if h5f.root.status._v_attrs.genomesAnnotated:
            logging.getLogger().info("Reading pangenome annotations...")
            readAnnotation(pangenome, h5f, filename)
            logging.getLogger().info("Done reading pangenome annotations")
        else:
            raise Exception(f"The pangenome in file '{filename}' has not been annotated, or has been improperly filled")
    if geneFamilies:
        if h5f.root.status._v_attrs.genesClustered:
            logging.getLogger().info("Reading pangenome gene families...")
            readGeneFamilies(pangenome, h5f)
            logging.getLogger().info("Done reading the pangenome's gene families")
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have gene families, or has been improperly filled")
    if graph:
        if h5f.root.status._v_attrs.NeighborsGraph:
            readGraph(pangenome, h5f)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have graph informations, or has been improperly filled")
    h5f.close()