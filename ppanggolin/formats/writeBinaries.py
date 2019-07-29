#!/usr/bin/env python3
#coding:utf-8

#default libraries
from pathlib import Path
import logging
import os
import warnings

#installed libraries
from tqdm import tqdm
import tables

#local libraries
from ppanggolin.pangenome import Pangenome, GeneFamily, Edge
from ppanggolin.genome import Organism, Contig, Gene
from ppanggolin.utils import getCurrentRAM

def geneDesc():
    return {
            'organism':tables.StringCol(itemsize=1024),
            "contig":{
                    'name':tables.StringCol(itemsize=1024),
                    "is_circular":tables.BoolCol()
            },
            "gene":{ 
                'ID':tables.StringCol(itemsize=1024),
                'start':tables.UInt64Col(),
                'stop':tables.UInt64Col(),
                'strand':tables.StringCol(itemsize=1),
                'type':tables.StringCol(itemsize=32),
                'position':tables.UInt32Col(),
                'name':tables.StringCol(itemsize=1024),
                'product':tables.StringCol(itemsize=1024),
                'genetic_code':tables.UInt32Col(),
                'is_fragment':tables.BoolCol(dflt = False),
            }
            }


def geneSequenceDesc():
    return {
        "gene":tables.StringCol(itemsize=1024),
        "dna":tables.StringCol(itemsize=160000),
        "type":tables.StringCol(itemsize=32)
    }

def geneFamDesc():
     return {
        "name": tables.StringCol(itemsize = 1024),
        "protein": tables.StringCol(itemsize=40000)
        }

def gene2famDesc():
    return {
        "geneFam": tables.StringCol(itemsize = 1024),
        "gene":tables.StringCol(itemsize= 1024)
        }

def writeAnnotations(pangenome, h5f):
    """
        Function writing all of the pangenome's annotations
    """
    # original_warnings = list(warnings.filters)
    # warnings.simplefilter('ignore', tables.NaturalNameWarning)
    annotation = h5f.create_group("/","annotations","Annotations of the pangenome's organisms")
    table = h5f.create_table(annotation, "genes", geneDesc(), expectedrows=len(pangenome.genes))
    bar = tqdm(pangenome.organisms, unit="genome")
    geneRow = table.row
    for org in bar:
        for contig in org.contigs:
            for gene in contig.genes:
                geneRow["organism"] = org.name
                geneRow["contig/name"] = contig.name
                geneRow["contig/is_circular"] = contig.is_circular#this should be somewhere else.
                geneRow["gene/ID"]= gene.ID
                geneRow["gene/start"] = gene.start
                geneRow["gene/stop"] = gene.stop
                geneRow["gene/strand"] = gene.strand
                geneRow["gene/type"] = gene.type
                geneRow["gene/position"] = gene.position
                geneRow["gene/name"] = gene.name
                geneRow["gene/product"] = gene.product
                geneRow["gene/is_fragment"] = gene.is_fragment
                geneRow["gene/genetic_code"] = gene.genetic_code
                geneRow.append()
    table.flush()

    bar.close()
    # warnings.filters = original_warnings

def writeGeneSequences(pangenome, h5f):
    geneSeq = h5f.create_table("/","geneSequences", geneSequenceDesc(), expectedrows=len(pangenome.genes))
    geneRow = geneSeq.row
    bar = tqdm(pangenome.genes, unit = "gene")
    for gene in bar:
        geneRow["gene"] = gene.ID
        geneRow["dna"] = gene.dna
        geneRow["type"] = gene.type
        geneRow.append()
    geneSeq.flush()
    bar.close()

def writeGeneFamSeq(pangenome, h5f):
    geneFamSeq = h5f.create_table("/","geneFamiliesSequences",geneFamDesc(), expectedrows=len(pangenome.geneFamilies))
    row = geneFamSeq.row
    bar = tqdm( pangenome.geneFamilies, unit = "gene family")
    for fam in bar:
        row["name"] = fam.name
        row["protein"] = fam.sequence
        row.append()
    geneFamSeq.flush()
    bar.close()

def writeGeneFamilies(pangenome, h5f):
    """
        Function writing all of the pangenome's gene families
    """
    geneFamilies = h5f.create_table("/", "geneFamilies",gene2famDesc())
    geneRow = geneFamilies.row
    bar = tqdm(pangenome.geneFamilies, unit = "gene family")
    for geneFam in bar:
        for gene in geneFam.genes:
            geneRow["gene"] = gene.ID
            geneRow["geneFam"] = geneFam.name
            geneRow.append()
    geneFamilies.flush()
    bar.close()

def writeGraph(pangenome, h5f):
    edges = h5f.create_group("/", "graph","Edges of the neighbors graph")
    counter = 0#for naming uniquely edges, as they don't really have a proper way of identifying themselves.
    for edge in pangenome.edges:
        edgeGroup = h5f.create_group(edges, str(counter))
        edgeGroup.source = edge.source.name
        edgeGroup.target = edge.target.name
        edgeTable = h5f.create_table(edgeGroup, 'genes', { 'organism':tables.StringCol(), 'geneTarget':tables.StringCol(), 'geneSource':tables.StringCol() })
        counter +=1
        edgeRow = edgeTable.row
        for org, genePairs in edge.organisms.items():
            for gene1, gene2 in genePairs:
                edgeRow["organism"] = org.name
                edgeRow["geneTarget"] = gene1.name
                edgeRow["geneSource"] = gene2.name
                edgeRow.append()
        edgeTable.flush()

def writeStatus(pangenome, h5f):
    if "/status" in h5f:#if statuses are already written
        statusGroup = h5f.root.status
    else:#else create the status group.
        statusGroup = h5f.create_group("/","status","Statuses of the pangenome's content")
    statusGroup._v_attrs.genomesAnnotated = True if pangenome.status["genomesAnnotated"] in ["Computed","Loaded","inFile"] else False
    statusGroup._v_attrs.geneSequences = True if pangenome.status["geneSequences"] in ["Computed","Loaded","inFile"] else False
    statusGroup._v_attrs.genesClustered = True if pangenome.status["genesClustered"] in ["Computed","Loaded","inFile"] else False
    statusGroup._v_attrs.geneFamilySequences = True if pangenome.status["geneFamilySequences"] in ["Computed","Loaded","inFile"] else False
    statusGroup._v_attrs.NeighborsGraph = True if pangenome.status["NeighborsGraph"] in ["Computed","Loaded","inFile"] else False
    statusGroup._v_attrs.Partitionned = True if pangenome.status["Partitionned"] in ["Computed","Loaded","inFile"] else False


def writePangenome(pangenome, filename):
    """
        Writes or updates a pangenome file
        pangenome is the corresponding pangenome object, filename the h5 file and status what has been modified.
    """

    compressionFilter = tables.Filters(complevel=1, complib='blosc:lz4')
    if pangenome.status["genomesAnnotated"] == "Computed":
        if filename.suffix != ".h5":
            filename = filename.with_suffix(".h5")
        h5f = tables.open_file(filename,"w", filters=compressionFilter)
        logging.getLogger().info("Writing genome annotations...")
        writeAnnotations(pangenome, h5f)
        h5f.close()
    elif pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"]:
        pass
        
    else:#if the pangenome is not Yes not Loaded, it's probably not really in a good state ( or something new was coded).
        raise NotImplementedError("Something REALLY unexpected and unplanned for happened here. Dev's contact is ggautrea [at] genoscope [dot] cns [dot] fr.")
    
    #from there, appending to existing file.
    h5f = tables.open_file(filename,"a", filters=compressionFilter)

    if pangenome.status["geneSequences"] == "Computed":
        logging.getLogger().info("writing the genes dna sequences")
        writeGeneSequences(pangenome, h5f)

    if pangenome.status["geneFamilySequences"] == "Computed":
        logging.getLogger().info("Writing the gene family representative protein sequences...")
        writeGeneFamSeq(pangenome, h5f)

    if pangenome.status["genesClustered"] == "Computed":
        logging.getLogger().info("Writing gene families and gene associations...")
        writeGeneFamilies(pangenome, h5f)
    
    if pangenome.status["NeighborsGraph"] == "Computed":
        raise NotImplementedError()
        # writeGraph(pangenome, h5f)
    
    if pangenome.status["Partitionned"] == "Computed":
        raise NotImplementedError()
        ##update geneFamilies with their partition.

    writeStatus(pangenome, h5f)
    h5f.close()
    logging.getLogger().info(f"Done writing the pangenome. It is in file : {filename}")