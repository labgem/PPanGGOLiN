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

def geneDesc(orgLen, contigLen, IDLen, typeLen, nameLen, productLen):
    return {
            'organism':tables.StringCol(itemsize=orgLen),
            "contig":{
                    'name':tables.StringCol(itemsize=contigLen),
                    "is_circular":tables.BoolCol(dflt = False)
            },
            "gene":{ 
                'ID':tables.StringCol(itemsize=IDLen),
                'start':tables.UInt64Col(),
                'stop':tables.UInt64Col(),
                'strand':tables.StringCol(itemsize=1),
                'type':tables.StringCol(itemsize=typeLen),
                'position':tables.UInt32Col(),
                'name':tables.StringCol(itemsize=nameLen),
                'product':tables.StringCol(itemsize=productLen),
                'genetic_code':tables.UInt32Col(dflt = 11),
                'is_fragment':tables.BoolCol(dflt = False),
            }
            }

def getMaxLenAnnotations(pangenome):
    maxOrgLen = 1
    maxContigLen = 1
    maxGeneIDLen = 1
    maxNameLen = 1
    maxProductLen = 1
    maxTypeLen = 1
    for org in pangenome.organisms:
        if len(org.name) > maxOrgLen:
            maxOrgLen = len(org.name)
        for contig in org.contigs:
            if len(contig.name) > maxContigLen:
                maxContigLen = len(contig.name)
            for gene in contig.genes:
                if len(gene.ID) > maxGeneIDLen:
                    maxGeneIDLen = len(gene.ID)
                if len(gene.name) > maxNameLen:
                    maxNameLen = len(gene.name)
                if len(gene.product) > maxProductLen:
                    maxProductLen = len(gene.product)
                if len(gene.type) > maxTypeLen:
                    maxTypeLen = len(gene.type)
    return maxOrgLen, maxContigLen, maxGeneIDLen, maxTypeLen, maxNameLen, maxProductLen

def writeAnnotations(pangenome, h5f):
    """
        Function writing all of the pangenome's annotations
    """
    annotation = h5f.create_group("/","annotations","Annotations of the pangenome's organisms")
    geneTable = h5f.create_table(annotation, "genes", geneDesc(*getMaxLenAnnotations(pangenome)), expectedrows=len(pangenome.genes))
    nbRNA = 0
    for org in pangenome.organisms:
        for contig in org.contigs:
            nbRNA += len(contig.RNAs)
    rnaTable = h5f.create_table(annotation, "RNA",geneDesc(*getMaxLenAnnotations(pangenome)), expectedrows=nbRNA)
    rnaRow = rnaTable.row
    bar = tqdm(pangenome.organisms, unit="genome")
    geneRow = geneTable.row
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
            for rna in contig.RNAs:
                rnaRow["organism"] = org.name
                rnaRow["contig/name"] = contig.name
                rnaRow["contig/is_circular"] = contig.is_circular#this should be somewhere else.
                rnaRow["gene/ID"]= rna.ID
                rnaRow["gene/start"] = rna.start
                rnaRow["gene/stop"] = rna.stop
                rnaRow["gene/strand"] = rna.strand
                rnaRow["gene/type"] = rna.type
                rnaRow["gene/name"] = rna.name
                rnaRow["gene/product"] = rna.product
                rnaRow["gene/is_fragment"] = rna.is_fragment
    geneTable.flush()
    rnaTable.flush()
    bar.close()


def getGeneSequencesLen(pangenome):
    maxSeqLen = 1
    maxGeneIDLen = 1
    maxGeneType = 1
    for gene in pangenome.genes:
        if len(gene.dna) > maxSeqLen:
            maxSeqLen = len(gene.dna)
        if len(gene.ID) > maxGeneIDLen:
            maxGeneIDLen = len(gene.ID)
        if len(gene.type) > maxGeneType:
            maxGeneType = len(gene.type)
    return maxGeneIDLen, maxSeqLen, maxGeneType

def geneSequenceDesc(geneIDLen, geneSeqLen, geneTypeLen):
    return {
        "gene":tables.StringCol(itemsize=geneIDLen),
        "dna":tables.StringCol(itemsize=geneSeqLen),
        "type":tables.StringCol(itemsize=geneTypeLen)
    }

def writeGeneSequences(pangenome, h5f):
    geneSeq = h5f.create_table("/","geneSequences", geneSequenceDesc(*getGeneSequencesLen(pangenome)), expectedrows=len(pangenome.genes))
    geneRow = geneSeq.row
    bar = tqdm(pangenome.genes, unit = "gene")
    for gene in bar:
        geneRow["gene"] = gene.ID
        geneRow["dna"] = gene.dna
        geneRow["type"] = gene.type
        geneRow.append()
    geneSeq.flush()
    bar.close()


def geneFamDesc(maxNameLen, maxSequenceLength):
     return {
        "name": tables.StringCol(itemsize = maxNameLen),
        "protein": tables.StringCol(itemsize=maxSequenceLength)
        }

def getGeneFamLen(pangenome):
    maxGeneFamNameLen = 1
    maxGeneFamSeqLen = 1
    for genefam in pangenome.geneFamilies:
        if len(genefam.sequence) > maxGeneFamSeqLen:
            maxGeneFamSeqLen = len(genefam.sequence)
        if len(genefam.name) > maxGeneFamNameLen:
            maxGeneFamNameLen = len(genefam.name)
    return maxGeneFamNameLen, maxGeneFamSeqLen

def writeGeneFamSeq(pangenome, h5f, force):
    """
        Writing a table containing the protein sequences of each family
    """
    if '/geneFamiliesSequences' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed gene family representative sequences...")
        h5f.remove_node('/', 'geneFamiliesSequences')#erasing the table, and rewriting a new one.
    geneFamSeq = h5f.create_table("/","geneFamiliesSequences",geneFamDesc(*getGeneFamLen(pangenome)), expectedrows=len(pangenome.geneFamilies))
    row = geneFamSeq.row
    bar = tqdm( pangenome.geneFamilies, unit = "gene family")
    for fam in bar:
        row["name"] = fam.name
        row["protein"] = fam.sequence
        row.append()
    geneFamSeq.flush()
    bar.close()


def gene2famDesc(geneFamNameLen, geneIDLen):
    return {
        "geneFam": tables.StringCol(itemsize = geneFamNameLen),
        "gene":tables.StringCol(itemsize= geneIDLen)
        }

def getGene2famLen(pangenome):
    maxGeneFamName = 1
    maxGeneID = 1
    for geneFam in pangenome.geneFamilies:
        if len(geneFam.name)>maxGeneFamName:
            maxGeneFamName = len(geneFam.name)
        for gene in geneFam.genes:
            if len(gene.ID) > maxGeneID:
                maxGeneID = len(gene.ID)
    return maxGeneFamName, maxGeneID

def writeGeneFamilies(pangenome, h5f, force):
    """
        Function writing all of the pangenome's gene families
    """
    if '/geneFamilies' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed gene family to gene associations...")
        h5f.remove_node('/', 'geneFamilies')#erasing the table, and rewriting a new one.
    geneFamilies = h5f.create_table("/", "geneFamilies",gene2famDesc(*getGene2famLen(pangenome)))
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
    statusGroup._v_attrs.defragmented = True if pangenome.status["defragmented"] in ["Computed","Loaded","inFile"] else False


def updateGeneFragments(pangenome, h5f):
    """
        updates the annotation table with the fragmentation informations from the defrag pipeline
    """
    logging.getLogger().info("Updating annotations with fragment informations")
    table = h5f.root.annotations.genes
    row = table.row
    bar =  tqdm(range(table.nrows), unit="gene")
    for row in table:
        row['gene/is_fragment'] = pangenome.getGene(row['gene/ID'].decode()).is_fragment
        bar.update()
    bar.close()


def writePangenome(pangenome, filename, force):
    """
        Writes or updates a pangenome file
        pangenome is the corresponding pangenome object, filename the h5 file and status what has been modified.
    """

    compressionFilter = tables.Filters(complevel=1, complib='blosc:lz4')#test the other compressors from blosc, this one was arbitrarily chosen.
    if pangenome.status["genomesAnnotated"] == "Computed":
        h5f = tables.open_file(filename,"w", filters=compressionFilter)
        logging.getLogger().info("Writing genome annotations...")
        writeAnnotations(pangenome, h5f)
        h5f.close()
    elif pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"]:
        pass
        
    else:#if the pangenome is not Computed not Loaded, it's probably not really in a good state ( or something new was coded).
        raise NotImplementedError("Something REALLY unexpected and unplanned for happened here. Dev's contact is ggautrea [at] genoscope [dot] cns [dot] fr.")
    
    #from there, appending to existing file.
    h5f = tables.open_file(filename,"a", filters=compressionFilter)

    if pangenome.status["geneSequences"] == "Computed":
        logging.getLogger().info("writing the protein coding genes dna sequences")
        writeGeneSequences(pangenome, h5f)

    if pangenome.status["geneFamilySequences"] == "Computed":
        logging.getLogger().info("Writing the gene family representative protein sequences...")
        writeGeneFamSeq(pangenome, h5f, force)

    if pangenome.status["genesClustered"] == "Computed":
        logging.getLogger().info("Writing gene families and gene associations...")
        writeGeneFamilies(pangenome, h5f, force)
        if pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"] and pangenome.status["defragmented"] == "Computed":
            #if the annotations have not been computed in this run, and there has been a clustering with defragmentation, then the annotations can be updated
            updateGeneFragments(pangenome,h5f)
    if pangenome.status["NeighborsGraph"] == "Computed":
        raise NotImplementedError()
        # writeGraph(pangenome, h5f)
    
    if pangenome.status["Partitionned"] == "Computed":
        raise NotImplementedError()
        ##update geneFamilies with their partition.

    writeStatus(pangenome, h5f)
    h5f.close()
    logging.getLogger().info(f"Done writing the pangenome. It is in file : {filename}")