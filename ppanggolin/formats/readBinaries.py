#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import sys

#installed libraries
from tqdm import tqdm
import tables

#local libraries
from ppanggolin.genome import Organism, Gene, RNA


def getNumberOfOrganisms(pangenome):
    """ standalone function to get the number of organisms in a pangenome"""
    if hasattr(pangenome,"file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f =  tables.open_file(filename,"r")
    annotations = h5f.root.annotations

    table = annotations.genes
    orgSet = set()
    for org in read_chunks(table, column = "organism"):
        orgSet.add(org)
    h5f.close()
    return len(orgSet)

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
        pangenome.status["neighborsGraph"] = "inFile"
    if statusGroup._v_attrs.Partitionned:
        pangenome.status["partitionned"] = "inFile"
    
    if hasattr(statusGroup._v_attrs, "predictedRGP"):
        if statusGroup._v_attrs.predictedRGP:
            pangenome.status["predictedRGP"] = "inFile"

    if "/info" in h5f:
        infoGroup = h5f.root.info
        pangenome.parameters = infoGroup._v_attrs.parameters
    h5f.close()

def read_chunks(table, column = None, chunk=10000):
    """
        Reading entirely the provided table (or column if specified) chunk per chunk to limit RAM usage.
    """
    for i in range(0,  table.nrows, chunk):
        for row in table.read(start = i, stop = i + chunk, field = column):
            yield row

def getGeneSequencesFromFile(pangenome, fileObj, list_CDS=None):
    """
        Writes the CDS sequences of the Pangenome object to a tmpFile object that can by filtered or not by a list of CDS
        Loads the sequences from a .h5 pangenome file
    """
    logging.getLogger().info("Extracting and writing all of the CDS sequences from a .h5 pangenome file to a fasta file")
    h5f = tables.open_file(pangenome.file,"r", driver_core_backing_store=0)
    table = h5f.root.geneSequences
    bar =  tqdm(range(table.nrows), unit="gene")
    list_CDS=set(list_CDS) if list_CDS is not None else None
    for row in read_chunks(table, chunk = 20000):#reading the table chunk per chunk otherwise RAM dies on big pangenomes
        nameCDS = row[1].decode()
        if row[2] == b"CDS" and (list_CDS is None or nameCDS in list_CDS):
            fileObj.write('>' + nameCDS + "\n")
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
        contig = org.getOrAddContig(contigName, is_circular=geneList[0][0][0])
        for row in geneList:
            if link:#if the gene families are already computed/loaded the gene exists.
                gene = pangenome.getGene(row[1][0].decode())
            else:#else creating the gene.
                if row[1][9].decode() == "CDS":
                    gene = Gene(row[1][0].decode())
                elif "RNA" in row[1][9].decode().upper():
                    gene = RNA(row[1][0].decode())
            gene.fill_parents(org, contig)
            if row[1][9].decode() == "CDS":
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
                contig.addGene(gene)
            elif "RNA" in row[1][9].decode():
                contig.addRNA(gene)
                gene.fill_annotations(
                        start = row[1][6],
                        stop =row[1][7],
                        strand =  row[1][8].decode(),
                        geneType = row[1][9].decode(),
                        name = row[1][3].decode(),
                        product = row[1][5].decode())
            else:
                raise Exception(f"A strange type ({row[1][9].decode()}), which we do not know what to do with, was met.")
    pangenome.addOrganism(org)

def readGraph(pangenome, h5f):
    table = h5f.root.edges

    if not pangenome.status["genomesAnnotated"] in ["Computed","Loaded"] or not pangenome.status["genesClustered"] in ["Computed","Loaded"] :
        raise Exception("It's not possible to read the graph if the annotations and the gene families have not been loaded.")
    bar = tqdm(range(table.nrows), unit = "contig adjacency")
    for row in read_chunks(table):
        source = pangenome.getGene(row[0].decode())
        target = pangenome.getGene(row[1].decode())
        pangenome.addEdge(source, target)
        bar.update()
    bar.close()
    pangenome.status["neighborsGraph"] = "Loaded"

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

def readGeneFamiliesInfo(pangenome, h5f):
    table = h5f.root.geneFamiliesInfo

    bar = tqdm(range(table.nrows), unit = "gene family")
    for row in read_chunks(table):
        fam = pangenome.addGeneFamily(row[0].decode())
        fam.addPartition(row[1].decode())
        fam.addSequence(row[2].decode())
        bar.update()
    bar.close()
    if h5f.root.status._v_attrs.Partitionned:
        pangenome.status["partitionned"] = "Loaded"
    if h5f.root.status._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "Loaded"


def readRGP(pangenome, h5f):
    table = h5f.root.RGP

    bar = tqdm(range(table.nrows), unit = "gene")
    for row in read_chunks(table):
        region = pangenome.getOrAddRegion(row[0].decode())
        region.append(pangenome.getGene(row[1].decode()))
        bar.update()
    bar.close()
    #order the genes properly in the regions
    for region in pangenome.regions:
        region.genes = sorted(region.genes, key = lambda x : x.position )#order the same way than on the contig
    pangenome.status["predictedRGP"] = "Loaded"

def readAnnotation(pangenome, h5f, filename):
    annotations = h5f.root.annotations

    table = annotations.genes
    bar = tqdm(range(table.nrows), unit="gene")
    pangenomeDict = {}
    for row in read_chunks(table):
        try:
            pangenomeDict[row[2].decode()][row[0][1].decode()].append(row)#new gene, seen contig, seen org
        except KeyError:
            try:
                pangenomeDict[row[2].decode()][row[0][1].decode()] = [row]#new contig, seen org
            except KeyError:
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

def readInfo(h5f):
    if "/info" in h5f:
        infoGroup = h5f.root.info

        print(f"Genes : {infoGroup._v_attrs['numberOfGenes']}")
        if "numberOfOrganisms" in infoGroup._v_attrs._f_list():
            print(f"Organisms : {infoGroup._v_attrs['numberOfOrganisms']}")
        if "numberOfClusters" in infoGroup._v_attrs._f_list():
            print(f"Families : {infoGroup._v_attrs['numberOfClusters']}")
        if "numberOfEdges" in infoGroup._v_attrs._f_list():
            print(f"Edges : {infoGroup._v_attrs['numberOfEdges']}")
        if 'numberOfCloud' in infoGroup._v_attrs._f_list():#then all the others are there
            print(f"Persistent ( { ', '.join([key + ':' + str(round(val,2)) for key, val in infoGroup._v_attrs['persistentStats'].items()])} ): {infoGroup._v_attrs['numberOfPersistent']}")
            print(f"Shell ( { ', '.join([key + ':' + str(round(val,2)) for key, val in infoGroup._v_attrs['shellStats'].items()])} ): {infoGroup._v_attrs['numberOfShell']}")
            print(f"Cloud ( { ', '.join([key + ':' + str(round(val,2)) for key, val in infoGroup._v_attrs['cloudStats'].items()])} ): {infoGroup._v_attrs['numberOfCloud']}")
            print(f"Number of partitions : {infoGroup._v_attrs['numberOfPartitions']}")
            if infoGroup._v_attrs['numberOfPartitions'] != 3:
                for key, val in infoGroup._v_attrs['numberOfSubpartitions'].items():
                    print(f"Shell {key} : {val}")
        if 'numberOfRGP' in infoGroup._v_attrs._f_list():
            print(f"RGPs : {infoGroup._v_attrs['numberOfRGP']}")

def readParameters(h5f):
    if "/info" in h5f:
        infoGroup = h5f.root.info
        if "parameters" in infoGroup._v_attrs._f_list():
            for key, dic in infoGroup._v_attrs["parameters"].items():
                print(f"{key}")
                for key2, val in dic.items():
                    print(f"    {key2} : {val}")

def readPangenome(pangenome, annotation = False, geneFamilies = False, graph = False, rgp = False):
    """
        Reads a previously written pangenome, with all of its parts, depending on what is asked.
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
        else:
            raise Exception(f"The pangenome in file '{filename}' has not been annotated, or has been improperly filled")
    if geneFamilies:
        if h5f.root.status._v_attrs.genesClustered:
            logging.getLogger().info("Reading pangenome gene families...")
            readGeneFamilies(pangenome, h5f)
            readGeneFamiliesInfo(pangenome, h5f)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have gene families, or has been improperly filled")
    if graph:
        if h5f.root.status._v_attrs.NeighborsGraph:
            logging.getLogger().info("Reading the neighbors graph edges...")
            readGraph(pangenome, h5f)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have graph informations, or has been improperly filled")
    if rgp:
        if h5f.root.status._v_attrs.predictedRGP:
            logging.getLogger().info("Reading the RGP...")
            readRGP(pangenome, h5f)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have RGP informations, or has been improperly filled")
    h5f.close()

def checkPangenomeInfo(pangenome, needAnnotations = False, needFamilies = False, needGraph = False, needPartitions = False, needRGP = False):
    """defines what needs to be read depending on what is needed"""
    annotation = False
    geneFamilies = False
    graph = False
    rgp = False
    if needAnnotations:
        if pangenome.status["genomesAnnotated"] == "inFile":
            annotation = True
        elif not pangenome.status["genomesAnnotated"] in ["Computed","Loaded"]:
            raise Exception("Your pangenome has no genes. See the 'annotate' subcommand.")
    if needFamilies:
        if pangenome.status["genesClustered"] == "inFile":
            geneFamilies = True
        elif not pangenome.status["genesClustered"] in ["Computed","Loaded"]:
            raise Exception("Your pangenome has no gene families. See the 'cluster' subcommand.")
    if needGraph:
        if pangenome.status["neighborsGraph"] == "inFile":
            graph=True
        elif not pangenome.status["neighborsGraph"] in ["Computed","Loaded"]:
            raise Exception("Your pangenome does not have a graph (no edges). See the 'graph' subcommand.")
    if needPartitions:
        if not pangenome.status["partitionned"] in ["Computed", "Loaded", "inFile"]:
            raise Exception("Your pangenome has not been partitionned. See the 'partition' subcommand")
    if needRGP:
        if pangenome.status["predictedRGP"] == "inFile":
            rgp = True
        elif not pangenome.status["predictedRGP"] in ["Computed","Loaded"]:
            raise Exception("Your pangenome's regions of genomic plasticity have not been predicted. See the 'rgp' subcommand")
    if annotation or geneFamilies or graph:#if anything is true, else we need nothing.
        readPangenome(pangenome, annotation=annotation, geneFamilies=geneFamilies, graph=graph, rgp = rgp)