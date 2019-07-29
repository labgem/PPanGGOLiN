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
from ppanggolin.annotate import genetic_codes, translate


def writeCDSSequences(pangenome, fileObj, cpu):
    logging.getLogger().info("Extracting and writing all of the CDS sequences")
    h5f = tables.open_file(pangenome,"r")
    table = h5f.root.annotations.genes
    bar =  tqdm(range(table.nrows), unit="gene")
    for row in table.read():
        if row[1][10].decode() == "CDS":
            fileObj.write('>' + row[1][0].decode() + "\n")
            fileObj.write(row[1][1].decode() + "\n")
        bar.update()

    bar.close()
    h5f.close()

def launchReadOrganism(args):
    return readOrganism(*args)

def readOrganism(orgName, contigDict, withSeq):
    org = Organism(orgName)
    for contigName, geneList in contigDict.items():
        contig = org.addContig(contigName, is_circular=geneList[0][0][0])
        for row in geneList:
            gene = Gene(row[1][0])
            if withSeq is True:
                gene.add_dna(row[1][1].decode())
            gene.fill_annotations(
                start = row[1][7],
                stop =row[1][8],
                strand =  row[1][9].decode(),
                geneType = row[1][10].decode(),
                position = row[1][5],
                genetic_code=row[1][2],
                name = row[1][4].decode(),
                product = row[1][6].decode())
            gene.is_fragment = row[1][3]
            if gene.type == "CDS":
                contig.addGene(gene)
            elif "RNA" in gene.type.upper():
                contig.addRNA(gene)
            else:
                raise Exception(f"A strange type ({gene.type}), which we do not know what to do with, was met.")
    return org

def readGraph(pangenome, h5f):
    raise NotImplementedError()

def readGeneFamilies(pangenome, h5f, partitions=True, genes = True, organisms = True):#reads everything ... but only if it is there.
    raise NotImplementedError()

def readAnnotation(pangenome, h5f, filename, cpu, withSeq = True):
    annotations = h5f.root.annotations
    # h5fpaths = []
    
    table = annotations.genes
    bar = tqdm(range(table.nrows), unit="gene")
    pangenomeDict = {}
    for row in table.read():
        try:
            pangenomeDict[row[2].decode()][row[0][1].decode()].append(row)#new gene, seen contig, seen org
        except:
            try:
                pangenomeDict[row[2].decode()][row[0][1].decode()] = [row]#new contig, seen org
            except:
                pangenomeDict[sys.intern(row[2].decode())] = { row[0][1].decode() : [row]}#new org
        bar.update()
    bar.close()
    bar = tqdm(range(len(pangenomeDict)), unit = "organism")
    for orgName, contigDict in pangenomeDict.items():
        pangenome.addOrganism(readOrganism(orgName, contigDict, withSeq))
        bar.update()
    bar.close()

def readPangenome(filename, cpu = 1, annotation = False, geneFamilies = False, graph = False):
    """
        Reads a previously written pangenome, with all of its parts.
    """
    pangenome = Pangenome()
    # compressionFilter = tables.Filters(complevel=1, complib='blosc:lz4')
    h5f = tables.open_file(filename,"r")
    if annotation:
        if h5f.root.status._v_attrs.genomesAnnotated:
            logging.getLogger().info("Reading pangenome annotations...")
            readAnnotation(pangenome, h5f, filename, cpu)
            logging.getLogger().info("Done reading pangenome annotations")
        else:
            raise Exception(f"The pangenome in file '{filename}' has not been annotated, or has been improperly filled")
    if geneFamilies:
        if h5f.root.status._v_attrs.genesClustered:
            readGeneFamilies(pangenome, h5f)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have gene families, or has been improperly filled")
    if graph:
        if h5f.root.status._v_attrs.NeighborsGraph:
            readGraph(pangenome, h5f)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have graph informations, or has been improperly filled")
    h5f.close()
    return pangenome