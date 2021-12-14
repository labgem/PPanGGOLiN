#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import sys

# installed libraries
from tqdm import tqdm
import tables

# local libraries
from ppanggolin.genome import Organism, Gene, RNA
from ppanggolin.region import Spot, Module


def getNumberOfOrganisms(pangenome):
    """ standalone function to get the number of organisms in a pangenome"""
    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")
    annotations = h5f.root.annotations

    table = annotations.genes
    orgSet = set()
    for org in read_chunks(table, column="organism"):
        orgSet.add(org)
    h5f.close()
    return len(orgSet)


def getStatus(pangenome, pangenomeFile):
    """
        Checks which elements are already present in the file.
    """
    h5f = tables.open_file(pangenomeFile, "r")
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

    if hasattr(statusGroup._v_attrs, "predictedRGP") and statusGroup._v_attrs.predictedRGP:
        pangenome.status["predictedRGP"] = "inFile"

    if hasattr(statusGroup._v_attrs, "spots") and statusGroup._v_attrs.spots:
        pangenome.status["spots"] = "inFile"

    if hasattr(statusGroup._v_attrs, "modules") and statusGroup._v_attrs.modules:
        pangenome.status["modules"] = "inFile"

    if "/info" in h5f:
        infoGroup = h5f.root.info
        pangenome.parameters = infoGroup._v_attrs.parameters
    h5f.close()


def read_chunks(table, column=None, chunk=10000):
    """
        Reading entirely the provided table (or column if specified) chunk per chunk to limit RAM usage.
    """
    for i in range(0, table.nrows, chunk):
        for row in table.read(start=i, stop=i + chunk, field=column):
            yield row


def getGeneSequencesFromFile(filename, fileObj, list_CDS=None, add='', disable_bar=False):
    """
        Writes the CDS sequences of the Pangenome object to a File object that can be filtered or not by a list of CDS, and adds the eventual str 'add' in front of the identifiers
        Loads the sequences from a .h5 pangenome file
    """
    logging.getLogger().info("Extracting and writing CDS sequences from a .h5 pangenome file to a fasta file...")
    h5f = tables.open_file(filename, "r", driver_core_backing_store=0)
    table = h5f.root.geneSequences
    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    list_CDS = set(list_CDS) if list_CDS is not None else None
    for row in read_chunks(table,
                           chunk=20000):  # reading the table chunk per chunk otherwise RAM dies on big pangenomes
        nameCDS = row["gene"].decode()
        if row["type"] == b"CDS" and (list_CDS is None or nameCDS in list_CDS):
            fileObj.write('>' + add + nameCDS + "\n")
            fileObj.write(row["dna"].decode() + "\n")
        bar.update()
    fileObj.flush()
    bar.close()
    h5f.close()


def launchReadOrganism(args):
    return readOrganism(*args)


def readOrganism(pangenome, orgName, contigDict, circularContigs, link=False):
    org = Organism(orgName)
    gene, gene_type = (None, None)
    for contigName, geneList in contigDict.items():
        contig = org.getOrAddContig(contigName, is_circular=circularContigs[contigName])
        for row in geneList:
            if link:  # if the gene families are already computed/loaded the gene exists.
                gene = pangenome.getGene(row["ID"].decode())
            else:  # else creating the gene.
                gene_type = row["type"].decode()
                if gene_type == "CDS":
                    gene = Gene(row["ID"].decode())
                elif "RNA" in gene_type:
                    gene = RNA(row["ID"].decode())
            try:
                local = row["local"].decode()
            except ValueError:
                local = ""
            gene.fill_annotations(
                start=row["start"],
                stop=row["stop"],
                strand=row["strand"].decode(),
                geneType=row["type"].decode(),
                position=row["position"],
                genetic_code=row["genetic_code"],
                name=row["name"].decode(),
                product=row["product"].decode(),
                local_identifier=local)
            gene.is_fragment = row["is_fragment"]
            gene.fill_parents(org, contig)
            if gene_type == "CDS":
                contig.addGene(gene)
            elif "RNA" in gene_type:
                contig.addRNA(gene)
            else:
                raise Exception(f"A strange type '{gene_type}', which we do not know what to do with, was met.")
    pangenome.addOrganism(org)


def readGraph(pangenome, h5f, disable_bar=False):
    table = h5f.root.edges

    if not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] or \
            not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the graph "
                        "if the annotations and the gene families have not been loaded.")

    bar = tqdm(range(table.nrows), unit="contig adjacency", disable=disable_bar)
    for row in read_chunks(table):
        source = pangenome.getGene(row["geneSource"].decode())
        target = pangenome.getGene(row["geneTarget"].decode())
        pangenome.addEdge(source, target)
        bar.update()
    bar.close()
    pangenome.status["neighborsGraph"] = "Loaded"


def readGeneFamilies(pangenome, h5f, disable_bar=False):
    table = h5f.root.geneFamilies

    link = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] else False

    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    for row in read_chunks(table):
        fam = pangenome.addGeneFamily(row["geneFam"].decode())
        if link:  # linking if we have loaded the annotations
            geneObj = pangenome.getGene(row["gene"].decode())
        else:  # else, no
            geneObj = Gene(row["gene"].decode())
        fam.addGene(geneObj)
        bar.update()
    bar.close()
    pangenome.status["genesClustered"] = "Loaded"


def readGeneFamiliesInfo(pangenome, h5f, disable_bar=False):
    table = h5f.root.geneFamiliesInfo

    bar = tqdm(range(table.nrows), unit="gene family", disable=disable_bar)
    for row in read_chunks(table):
        fam = pangenome.addGeneFamily(row["name"].decode())
        fam.addPartition(row["partition"].decode())
        fam.addSequence(row["protein"].decode())
        bar.update()
    bar.close()
    if h5f.root.status._v_attrs.Partitionned:
        pangenome.status["partitionned"] = "Loaded"
    if h5f.root.status._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "Loaded"


def readGeneSequences(pangenome, h5f, disable_bar=False):
    table = h5f.root.geneSequences

    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    for row in read_chunks(table):
        gene = pangenome.getGene(row['gene'].decode())
        gene.add_dna(row['dna'].decode())
        bar.update()
    bar.close()
    pangenome.status["geneSequences"] = "Loaded"


def readRGP(pangenome, h5f, disable_bar=False):
    table = h5f.root.RGP

    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    for row in read_chunks(table):
        region = pangenome.getOrAddRegion(row["RGP"].decode())
        region.append(pangenome.getGene(row["gene"].decode()))
        bar.update()
    bar.close()
    # order the genes properly in the regions
    for region in pangenome.regions:
        region.genes = sorted(region.genes, key=lambda x: x.position)  # order the same way as on the contig
    pangenome.status["predictedRGP"] = "Loaded"


def readSpots(pangenome, h5f, disable_bar=False):
    table = h5f.root.spots
    bar = tqdm(range(table.nrows), unit="region", disable=disable_bar)
    spots = {}
    for row in read_chunks(table):
        curr_spot = spots.get(row["spot"])
        if curr_spot is None:
            curr_spot = Spot(row["spot"])
            spots[row["spot"]] = curr_spot
        curr_spot.addRegion(pangenome.getOrAddRegion(row["RGP"].decode()))
        bar.update()
    bar.close()
    pangenome.addSpots(spots.values())
    pangenome.status["spots"] = "Loaded"


def readModules(pangenome, h5f, disable_bar=False):
    table = h5f.root.modules
    bar = tqdm(range(table.nrows), unit="module", disable=disable_bar)
    modules = {}  # id2mod
    for row in read_chunks(table):
        curr_module = modules.get(row['module'])
        if curr_module is None:
            curr_module = Module(row['module'])
            modules[row["module"]] = curr_module
        curr_module.addFamily(pangenome.getGeneFamily(row['geneFam'].decode()))
        bar.update()
    bar.close()
    pangenome.addModules(modules.values())
    pangenome.status["modules"] = "Loaded"


def readAnnotation(pangenome, h5f, disable_bar=False):
    annotations = h5f.root.annotations

    table = annotations.genes
    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    pangenomeDict = {}
    circularContigs = {}
    for row in read_chunks(table):
        try:
            pangenomeDict[row["organism"].decode()][row["contig"]["name"].decode()].append(
                row["gene"])  # new gene, seen contig, seen org
        except KeyError:
            try:
                pangenomeDict[row["organism"].decode()][row["contig"]["name"].decode()] = [
                    row["gene"]]  # new contig, seen org
                circularContigs[row["organism"].decode()][row["contig"]["name"].decode()] = row["contig"]["is_circular"]
            except KeyError:
                pangenomeDict[sys.intern(row["organism"].decode())] = {
                    row["contig"]["name"].decode(): [row["gene"]]}  # new org
                circularContigs[row["organism"].decode()] = {
                    row["contig"]["name"].decode(): row["contig"]["is_circular"]}
        bar.update()
    bar.close()

    link = True if pangenome.status["genesClustered"] in ["Computed", "Loaded"] else False

    bar = tqdm(range(len(pangenomeDict)), unit="organism", disable=disable_bar)
    for orgName, contigDict in pangenomeDict.items():
        readOrganism(pangenome, orgName, contigDict, circularContigs[orgName], link)
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
        if 'numberOfCloud' in infoGroup._v_attrs._f_list():  # then all the others are there
            print(
                f"Persistent ( {', '.join([key + ':' + str(round(val, 2)) for key, val in infoGroup._v_attrs['persistentStats'].items()])} ): {infoGroup._v_attrs['numberOfPersistent']}")
            print(
                f"Shell ( {', '.join([key + ':' + str(round(val, 2)) for key, val in infoGroup._v_attrs['shellStats'].items()])} ): {infoGroup._v_attrs['numberOfShell']}")
            print(
                f"Cloud ( {', '.join([key + ':' + str(round(val, 2)) for key, val in infoGroup._v_attrs['cloudStats'].items()])} ): {infoGroup._v_attrs['numberOfCloud']}")
            print(f"Number of partitions : {infoGroup._v_attrs['numberOfPartitions']}")
            if infoGroup._v_attrs['numberOfPartitions'] != 3:
                for key, val in infoGroup._v_attrs['numberOfSubpartitions'].items():
                    print(f"Shell {key} : {val}")
        if 'numberOfRGP' in infoGroup._v_attrs._f_list():
            print(f"RGPs : {infoGroup._v_attrs['numberOfRGP']}")
        if 'numberOfSpots' in infoGroup._v_attrs._f_list():
            print(f"Spots : {infoGroup._v_attrs['numberOfSpots']}")
        if 'numberOfModules' in infoGroup._v_attrs._f_list():
            print(f"Modules : {infoGroup._v_attrs['numberOfModules']}")
            print(f"Families in Modules : {infoGroup._v_attrs['numberOfFamiliesInModules']}")


def readParameters(h5f):
    if "/info" in h5f:
        infoGroup = h5f.root.info
        if "parameters" in infoGroup._v_attrs._f_list():
            for key, dic in infoGroup._v_attrs["parameters"].items():
                print(f"{key}")
                for key2, val in dic.items():
                    print(f"    {key2} : {val}")


def readPangenome(pangenome, annotation=False, geneFamilies=False, graph=False, rgp=False, spots=False,
                  geneSequences=False, modules=False, disable_bar=False):
    """
        Reads a previously written pangenome, with all of its parts, depending on what is asked,
        with regard to what is filled in the 'status' field of the hdf5 file.
    """
    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")
    if annotation:
        if h5f.root.status._v_attrs.genomesAnnotated:
            logging.getLogger().info("Reading pangenome annotations...")
            readAnnotation(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' has not been annotated, or has been improperly filled")
    if geneSequences:
        if h5f.root.status._v_attrs.geneSequences:
            logging.getLogger().info("Reading pangenome gene dna sequences...")
            readGeneSequences(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have gene sequences, "
                            f"or has been improperly filled")

    if geneFamilies:
        if h5f.root.status._v_attrs.genesClustered:
            logging.getLogger().info("Reading pangenome gene families...")
            readGeneFamilies(pangenome, h5f, disable_bar=disable_bar)
            readGeneFamiliesInfo(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have gene families, or has been improperly filled")
    if graph:
        if h5f.root.status._v_attrs.NeighborsGraph:
            logging.getLogger().info("Reading the neighbors graph edges...")
            readGraph(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have graph information, or has been improperly filled")
    if rgp:
        if h5f.root.status._v_attrs.predictedRGP:
            logging.getLogger().info("Reading the RGP...")
            readRGP(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have RGP information, or has been improperly filled")
    if spots:
        if h5f.root.status._v_attrs.spots:
            logging.getLogger().info("Reading the spots...")
            readSpots(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have spots information, or has been improperly filled")
    if modules:
        if h5f.root.status._v_attrs.modules:
            logging.getLogger().info("Reading the modules...")
            readModules(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have modules information, or has been improperly filled")
    h5f.close()


def checkPangenomeInfo(pangenome, needAnnotations=False, needFamilies=False, needGraph=False, needPartitions=False,
                       needRGP=False, needSpots=False, needGeneSequences=False, needModules=False, disable_bar=False):
    """
    defines what needs to be read depending on what is needed, and automatically checks if the required elements
    have been computed with regard to the pangenome.status
    """
    annotation = False
    geneFamilies = False
    graph = False
    rgp = False
    spots = False
    geneSequences = False
    modules = False

    if needAnnotations:
        if pangenome.status["genomesAnnotated"] == "inFile":
            annotation = True
        elif not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome has no genes. See the 'annotate' subcommand.")
    if needFamilies:
        if pangenome.status["genesClustered"] == "inFile":
            geneFamilies = True
        elif not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome has no gene families. See the 'cluster' subcommand.")
    if needGraph:
        if pangenome.status["neighborsGraph"] == "inFile":
            graph = True
        elif not pangenome.status["neighborsGraph"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome does not have a graph (no edges). See the 'graph' subcommand.")
    if needPartitions and pangenome.status["partitionned"] not in ["Computed", "Loaded", "inFile"]:
        raise Exception("Your pangenome has not been partitioned. See the 'partition' subcommand")
    if needRGP:
        if pangenome.status["predictedRGP"] == "inFile":
            rgp = True
        elif not pangenome.status["predictedRGP"] in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome's regions of genomic plasticity have not been predicted. See the 'rgp' subcommand")
    if needSpots:
        if pangenome.status["spots"] == "inFile":
            spots = True
        elif not pangenome.status["spots"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome spots of insertion have not been predicted. See the 'spot' subcommand")
    if needGeneSequences:
        if pangenome.status["geneSequences"] == "inFile":
            geneSequences = True
        elif not pangenome.status["geneSequences"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome does not include gene sequences. "
                            "This is possible only if you provided your own cluster file with the 'cluster' subcommand")

    if needModules:
        if pangenome.status["modules"] == "inFile":
            modules = True
        elif not pangenome.status["modules"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome modules have not been predicted. See the 'module' subcommand")

    if annotation or geneFamilies or graph or rgp or spots or geneSequences or modules:
        # if anything is true, else we need nothing.
        readPangenome(pangenome, annotation=annotation, geneFamilies=geneFamilies, graph=graph, rgp=rgp, spots=spots,
                      geneSequences=geneSequences, modules=modules, disable_bar=disable_bar)
