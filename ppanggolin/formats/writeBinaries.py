#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
from collections import Counter

#installed libraries
from tqdm import tqdm
import tables

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


def geneFamDesc(maxNameLen, maxSequenceLength, maxPartLen):
     return {
        "name": tables.StringCol(itemsize = maxNameLen),
        "protein": tables.StringCol(itemsize=maxSequenceLength),
        "partition": tables.StringCol(itemsize=maxPartLen)
        }

def getGeneFamLen(pangenome):
    maxGeneFamNameLen = 1
    maxGeneFamSeqLen = 1
    maxPartLen = 1
    for genefam in pangenome.geneFamilies:
        if len(genefam.sequence) > maxGeneFamSeqLen:
            maxGeneFamSeqLen = len(genefam.sequence)
        if len(genefam.name) > maxGeneFamNameLen:
            maxGeneFamNameLen = len(genefam.name)
        if len(genefam.partition) > maxPartLen:
            maxPartLen = len(genefam.partition)
    return maxGeneFamNameLen, maxGeneFamSeqLen, maxPartLen

def writeGeneFamInfo(pangenome, h5f, force):
    """
        Writing a table containing the protein sequences of each family
    """
    if '/geneFamiliesInfo' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed gene family representative sequences...")
        h5f.remove_node('/', 'geneFamiliesInfo')#erasing the table, and rewriting a new one.
    geneFamSeq = h5f.create_table("/","geneFamiliesInfo",geneFamDesc(*getGeneFamLen(pangenome)), expectedrows=len(pangenome.geneFamilies))

    row = geneFamSeq.row
    bar = tqdm( pangenome.geneFamilies, unit = "gene family")
    for fam in bar:
        row["name"] = fam.name
        row["protein"] = fam.sequence
        row["partition"] = fam.partition
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


def graphDesc(maxGeneIDLen):
    return {
            'geneTarget':tables.StringCol(itemsize = maxGeneIDLen),
            'geneSource':tables.StringCol(itemsize = maxGeneIDLen)
        }

def getGeneIDLen(pangenome):
    maxGeneLen = 1
    for gene in pangenome.genes:
        if len(gene.ID) > maxGeneLen:
            maxGeneLen=len(gene.ID)
    return maxGeneLen

def writeGraph(pangenome, h5f, force):
    #if we want to be able to read the graph without reading the annotations (because it is one of the most time consumming parts to read), it might be good to add the organism name in the table here.
    #for now, forcing the read of annotations.
    if '/edges' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed edges")
        h5f.remove_node("/","edges")
    edgeTable = h5f.create_table("/","edges", graphDesc(getGeneIDLen(pangenome)), expectedrows=len(pangenome.edges))
    edgeRow = edgeTable.row
    bar = tqdm(pangenome.edges, unit = "edge")
    for edge in bar:
        for _, genePairs in edge.organisms.items():
            for gene1, gene2 in genePairs:
                edgeRow["geneTarget"] = gene1.ID
                edgeRow["geneSource"] = gene2.ID
                edgeRow.append()
    bar.close()
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
    statusGroup._v_attrs.NeighborsGraph = True if pangenome.status["neighborsGraph"] in ["Computed","Loaded","inFile"] else False
    statusGroup._v_attrs.Partitionned = True if pangenome.status["partitionned"] in ["Computed","Loaded","inFile"] else False
    statusGroup._v_attrs.defragmented = True if pangenome.status["defragmented"] in ["Computed","Loaded","inFile"] else False


def writeInfo(pangenome, h5f):
    if "/info" in h5f:
        infoGroup = h5f.root.info
    else:
        infoGroup = h5f.create_group("/","info","Informations about the pangenome's content")
    if pangenome.status["genomesAnnotated"] in ["Computed","Loaded"]:
        infoGroup._v_attrs.numberOfGenes = len(pangenome.genes)
        infoGroup._v_attrs.numberOfOrganisms = len(pangenome.organisms)
    if pangenome.status["genesClustered"] in ["Computed","Loaded"]:
        infoGroup._v_attrs.numberOfClusters = len(pangenome.geneFamilies)
    if pangenome.status["neighborsGraph"] in ["Computed","Loaded"]:
        infoGroup._v_attrs.numberOfEdges = len(pangenome.edges)
    if pangenome.status["partitionned"] in ["Computed","Loaded"]:
        namedPartCounter = Counter()
        partSet = set()
        for fam in pangenome.geneFamilies:
            namedPartCounter[fam.namedPartition] +=1
            partSet.add(fam.partition)
        infoGroup._v_attrs.numberOfPersistent = namedPartCounter["persistent"]
        infoGroup._v_attrs.numberOfShell = namedPartCounter["shell"]
        infoGroup._v_attrs.numberOfCloud = namedPartCounter["cloud"]
        infoGroup._v_attrs.numberOfPartitions = len(partSet)

def updateGeneFamPartition(pangenome, h5f):
    logging.getLogger().info("Updating gene families with partition information")
    table = h5f.root.geneFamiliesInfo
    bar = tqdm(range(table.nrows), unit = "gene family")
    for row in table:
        row["partition"] = pangenome.getGeneFamily(row["name"].decode()).partition
        row.update()
        bar.update()
    bar.close()

def updateGeneFragments(pangenome, h5f):
    """
        updates the annotation table with the fragmentation informations from the defrag pipeline
    """
    logging.getLogger().info("Updating annotations with fragment information")
    table = h5f.root.annotations.genes
    row = table.row
    bar =  tqdm(range(table.nrows), unit="gene")
    for row in table:
        row['gene/is_fragment'] = pangenome.getGene(row['gene/ID'].decode()).is_fragment
        bar.update()
    bar.close()


def ErasePangenome(pangenome, graph=False, geneFamilies = False):
    """ erases tables from a pangenome .h5 file """
    compressionFilter = tables.Filters(complevel=1, complib='blosc:lz4')
    h5f = tables.open_file(pangenome.file,"a", filters=compressionFilter)
    statusGroup = h5f.root.status

    if '/edges' in h5f and (graph or geneFamilies):
        logging.getLogger().info("Erasing the formerly computed edges")
        h5f.remove_node("/","edges")
        statusGroup._v_attrs.NeighborsGraph = False
        pangenome.status["neighborsGraph"] = "No"
    if '/geneFamilies' in h5f and geneFamilies:
        logging.getLogger().info("Erasing the formerly computed gene family to gene associations...")
        h5f.remove_node('/', 'geneFamilies')#erasing the table, and rewriting a new one.
        pangenome.status["defragmented"] = "No"
        pangenome.status["genesClustered"] = "No"
        statusGroup._v_attrs.defragmented = False
        statusGroup._v_attrs.genesClustered = False
    if '/geneFamiliesInfo' in h5f and geneFamilies:
        logging.getLogger().info("Erasing the formerly computed gene family representative sequences...")
        h5f.remove_node('/', 'geneFamiliesInfo')#erasing the table, and rewriting a new one.
        pangenome.status["partitionned"] = "No"
        pangenome.status["geneFamilySequences"] = "No"
        statusGroup._v_attrs.geneFamilySequences = False
        statusGroup._v_attrs.Partitionned = False
    h5f.close()

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
        raise NotImplementedError("Something REALLY unexpected and unplanned for happened here. Please post an issue on github with what you did to reach this error.")

    #from there, appending to existing file.
    h5f = tables.open_file(filename,"a", filters=compressionFilter)

    if pangenome.status["geneSequences"] == "Computed":
        logging.getLogger().info("writing the protein coding gene dna sequences")
        writeGeneSequences(pangenome, h5f)

    if pangenome.status["genesClustered"] == "Computed":
        logging.getLogger().info("Writing gene families and gene associations...")
        writeGeneFamilies(pangenome, h5f, force)
        logging.getLogger().info("Writing gene families information...")
        writeGeneFamInfo(pangenome, h5f, force)
        if pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"] and pangenome.status["defragmented"] == "Computed":
            #if the annotations have not been computed in this run, and there has been a clustering with defragmentation, then the annotations can be updated
            updateGeneFragments(pangenome,h5f)
    if pangenome.status["neighborsGraph"] == "Computed":
        logging.getLogger().info("Writing the edges...")
        writeGraph(pangenome, h5f, force)

    if pangenome.status["partitionned"] == "Computed" and pangenome.status["genesClustered"] in ["Loaded","inFile"]:#otherwise it's been written already.
        updateGeneFamPartition(pangenome, h5f)

    writeStatus(pangenome, h5f)
    writeInfo(pangenome, h5f)

    h5f.close()
    logging.getLogger().info(f"Done writing the pangenome. It is in file : {filename}")