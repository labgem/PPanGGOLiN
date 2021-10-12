#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from collections import Counter, defaultdict
import statistics
import pkg_resources

# installed libraries
from tqdm import tqdm
import tables


def geneDesc(orgLen, contigLen, IDLen, typeLen, nameLen, productLen, maxLocalId):
    return {
        'organism': tables.StringCol(itemsize=orgLen),
        "contig": {
            'name': tables.StringCol(itemsize=contigLen),
            "is_circular": tables.BoolCol(dflt=False)
        },
        "gene": {
            'ID': tables.StringCol(itemsize=IDLen),
            'start': tables.UInt64Col(),
            'stop': tables.UInt64Col(),
            'strand': tables.StringCol(itemsize=1),
            'type': tables.StringCol(itemsize=typeLen),
            'position': tables.UInt32Col(),
            'name': tables.StringCol(itemsize=nameLen),
            'product': tables.StringCol(itemsize=productLen),
            'genetic_code': tables.UInt32Col(dflt=11),
            'is_fragment': tables.BoolCol(dflt=False),
            'local': tables.StringCol(itemsize=maxLocalId)
        }
    }


def getMaxLenAnnotations(pangenome):
    maxOrgLen = 1
    maxContigLen = 1
    maxGeneIDLen = 1
    maxNameLen = 1
    maxProductLen = 1
    maxTypeLen = 1
    maxLocalId = 1
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
                if len(gene.local_identifier) > maxLocalId:
                    maxLocalId = len(gene.local_identifier)
            for gene in contig.RNAs:
                if len(gene.ID) > maxGeneIDLen:
                    maxGeneIDLen = len(gene.ID)
                if len(gene.name) > maxNameLen:
                    maxNameLen = len(gene.name)
                if len(gene.product) > maxProductLen:
                    maxProductLen = len(gene.product)
                if len(gene.type) > maxTypeLen:
                    maxTypeLen = len(gene.type)
                if len(gene.local_identifier) > maxLocalId:
                    maxLocalId = len(gene.local_identifier)

    return maxOrgLen, maxContigLen, maxGeneIDLen, maxTypeLen, maxNameLen, maxProductLen, maxLocalId


def writeAnnotations(pangenome, h5f, disable_bar=False):
    """
        Function writing all of the pangenome's annotations
    """
    annotation = h5f.create_group("/", "annotations", "Annotations of the pangenome's organisms")
    geneTable = h5f.create_table(annotation, "genes", geneDesc(*getMaxLenAnnotations(pangenome)),
                                 expectedrows=len(pangenome.genes))
    nbRNA = 0
    for org in pangenome.organisms:
        for contig in org.contigs:
            nbRNA += len(contig.RNAs)
    bar = tqdm(pangenome.organisms, unit="genome", disable=disable_bar)
    geneRow = geneTable.row
    for org in bar:
        for contig in org.contigs:
            for gene in contig.genes:
                geneRow["organism"] = org.name
                geneRow["contig/name"] = contig.name
                geneRow["contig/is_circular"] = contig.is_circular  # this should be somewhere else.
                geneRow["gene/ID"] = gene.ID
                geneRow["gene/start"] = gene.start
                geneRow["gene/stop"] = gene.stop
                geneRow["gene/strand"] = gene.strand
                geneRow["gene/type"] = gene.type
                geneRow["gene/position"] = gene.position
                geneRow["gene/name"] = gene.name
                geneRow["gene/product"] = gene.product
                geneRow["gene/is_fragment"] = gene.is_fragment
                geneRow["gene/genetic_code"] = gene.genetic_code
                geneRow["gene/local"] = gene.local_identifier
                geneRow.append()
            for rna in contig.RNAs:
                geneRow["organism"] = org.name
                geneRow["contig/name"] = contig.name
                geneRow["contig/is_circular"] = contig.is_circular  # this should be somewhere else.
                geneRow["gene/ID"] = rna.ID
                geneRow["gene/start"] = rna.start
                geneRow["gene/stop"] = rna.stop
                geneRow["gene/strand"] = rna.strand
                geneRow["gene/type"] = rna.type
                geneRow["gene/name"] = rna.name
                geneRow["gene/product"] = rna.product
                geneRow["gene/is_fragment"] = rna.is_fragment
                geneRow.append()
    geneTable.flush()
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
        "gene": tables.StringCol(itemsize=geneIDLen),
        "dna": tables.StringCol(itemsize=geneSeqLen),
        "type": tables.StringCol(itemsize=geneTypeLen)
    }


def writeGeneSequences(pangenome, h5f, disable_bar=False):
    geneSeq = h5f.create_table("/", "geneSequences", geneSequenceDesc(*getGeneSequencesLen(pangenome)),
                               expectedrows=len(pangenome.genes))
    geneRow = geneSeq.row
    bar = tqdm(pangenome.genes, unit="gene", disable=disable_bar)
    for gene in bar:
        geneRow["gene"] = gene.ID
        geneRow["dna"] = gene.dna
        geneRow["type"] = gene.type
        geneRow.append()
    geneSeq.flush()
    bar.close()


def geneFamDesc(maxNameLen, maxSequenceLength, maxPartLen):
    return {
        "name": tables.StringCol(itemsize=maxNameLen),
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


def writeGeneFamInfo(pangenome, h5f, force, disable_bar=False):
    """
        Writing a table containing the protein sequences of each family
    """
    if '/geneFamiliesInfo' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed gene family representative sequences...")
        h5f.remove_node('/', 'geneFamiliesInfo')  # erasing the table, and rewriting a new one.
    geneFamSeq = h5f.create_table("/", "geneFamiliesInfo", geneFamDesc(*getGeneFamLen(pangenome)),
                                  expectedrows=len(pangenome.geneFamilies))

    row = geneFamSeq.row
    bar = tqdm(pangenome.geneFamilies, unit="gene family", disable=disable_bar)
    for fam in bar:
        row["name"] = fam.name
        row["protein"] = fam.sequence
        row["partition"] = fam.partition
        row.append()
    geneFamSeq.flush()
    bar.close()


def gene2famDesc(geneFamNameLen, geneIDLen):
    return {
        "geneFam": tables.StringCol(itemsize=geneFamNameLen),
        "gene": tables.StringCol(itemsize=geneIDLen)
    }


def getGene2famLen(pangenome):
    maxGeneFamName = 1
    maxGeneID = 1
    for geneFam in pangenome.geneFamilies:
        if len(geneFam.name) > maxGeneFamName:
            maxGeneFamName = len(geneFam.name)
        for gene in geneFam.genes:
            if len(gene.ID) > maxGeneID:
                maxGeneID = len(gene.ID)
    return maxGeneFamName, maxGeneID


def writeGeneFamilies(pangenome, h5f, force, disable_bar=False):
    """
        Function writing all of the pangenome's gene families
    """
    if '/geneFamilies' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed gene family to gene associations...")
        h5f.remove_node('/', 'geneFamilies')  # erasing the table, and rewriting a new one.
    geneFamilies = h5f.create_table("/", "geneFamilies", gene2famDesc(*getGene2famLen(pangenome)))
    geneRow = geneFamilies.row
    bar = tqdm(pangenome.geneFamilies, unit="gene family", disable=disable_bar)
    for geneFam in bar:
        for gene in geneFam.genes:
            geneRow["gene"] = gene.ID
            geneRow["geneFam"] = geneFam.name
            geneRow.append()
    geneFamilies.flush()
    bar.close()


def graphDesc(maxGeneIDLen):
    return {
        'geneTarget': tables.StringCol(itemsize=maxGeneIDLen),
        'geneSource': tables.StringCol(itemsize=maxGeneIDLen)
    }


def getGeneIDLen(pangenome):
    maxGeneLen = 1
    for gene in pangenome.genes:
        if len(gene.ID) > maxGeneLen:
            maxGeneLen = len(gene.ID)
    return maxGeneLen


def writeGraph(pangenome, h5f, force, disable_bar=False):
    # if we want to be able to read the graph without reading the annotations (because it is one of the most time consumming parts to read), it might be good to add the organism name in the table here.
    # for now, forcing the read of annotations.
    if '/edges' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed edges")
        h5f.remove_node("/", "edges")
    edgeTable = h5f.create_table("/", "edges", graphDesc(getGeneIDLen(pangenome)), expectedrows=len(pangenome.edges))
    edgeRow = edgeTable.row
    bar = tqdm(pangenome.edges, unit="edge", disable=disable_bar)
    for edge in bar:
        for genePairs in edge.organisms.values():
            for gene1, gene2 in genePairs:
                edgeRow["geneTarget"] = gene1.ID
                edgeRow["geneSource"] = gene2.ID
                edgeRow.append()
    bar.close()
    edgeTable.flush()


def RGPDesc(maxRGPLen, maxGeneLen):
    return {
        'RGP': tables.StringCol(itemsize=maxRGPLen),
        'gene': tables.StringCol(itemsize=maxGeneLen)
    }


def getRGPLen(pangenome):
    maxGeneLen = 1
    maxRGPLen = 1
    for region in pangenome.regions:
        for gene in region.genes:
            if len(gene.ID) > maxGeneLen:
                maxGeneLen = len(gene.ID)
        if len(region.name) > maxRGPLen:
            maxRGPLen = len(region.name)
    return maxRGPLen, maxGeneLen


def writeRGP(pangenome, h5f, force, disable_bar=False):
    if '/RGP' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computer RGP")
        h5f.remove_node('/', 'RGP')

    RGPTable = h5f.create_table('/', 'RGP', RGPDesc(*getRGPLen(pangenome)),
                                expectedrows=sum([len(region.genes) for region in pangenome.regions]))
    RGPRow = RGPTable.row
    bar = tqdm(pangenome.regions, unit="region", disable=disable_bar)
    for region in bar:
        for gene in region.genes:
            RGPRow["RGP"] = region.name
            RGPRow["gene"] = gene.ID
            RGPRow.append()
    bar.close()
    RGPTable.flush()


def spotDesc(maxRGPLen):
    return {
        'spot': tables.UInt32Col(),
        'RGP': tables.StringCol(itemsize=maxRGPLen)
    }


def getSpotDesc(pangenome):
    maxRGPLen = 1
    for spot in pangenome.spots:
        for region in spot.regions:
            if len(region.name) > maxRGPLen:
                maxRGPLen = len(region.name)
    return maxRGPLen


def writeSpots(pangenome, h5f, force, disable_bar=False):
    if '/spots' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed spots")
        h5f.remove_node("/", "spots")

    SpoTable = h5f.create_table("/", "spots", spotDesc(getSpotDesc(pangenome)),
                                expectedrows=sum([len(spot.regions) for spot in pangenome.spots]))
    SpotRow = SpoTable.row
    bar = tqdm(pangenome.spots, unit="spot", disable=disable_bar)
    for spot in pangenome.spots:
        for region in spot.regions:
            SpotRow["spot"] = spot.ID
            SpotRow["RGP"] = region.name
            SpotRow.append()
        bar.update()
    bar.close()
    SpoTable.flush()


def modDesc(geneFamNameLen):
    return {
        "geneFam": tables.StringCol(itemsize=geneFamNameLen),
        "module": tables.UInt32Col(),
    }


def getModDesc(pangenome):
    maxFamLen = 1
    for mod in pangenome.modules:
        for fam in mod.families:
            if len(fam.name) > maxFamLen:
                maxFamLen = len(fam.name)
    return maxFamLen


def writeModules(pangenome, h5f, force, disable_bar=False):
    if '/modules' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed modules")
        h5f.remove_node("/", "modules")

    modTable = h5f.create_table('/', 'modules', modDesc(getModDesc(pangenome)),
                                expectedrows=sum([len(mod.families) for mod in pangenome.modules]))
    modRow = modTable.row

    bar = tqdm(pangenome.modules, unit="modules", disable=disable_bar)
    for mod in bar:
        for fam in mod.families:
            modRow["geneFam"] = fam.name
            modRow["module"] = mod.ID
            modRow.append()
    bar.close()
    modTable.flush()


def writeStatus(pangenome, h5f):
    if "/status" in h5f:  # if statuses are already written
        statusGroup = h5f.root.status
    else:  # else create the status group.
        statusGroup = h5f.create_group("/", "status", "Statuses of the pangenome's content")
    statusGroup._v_attrs.genomesAnnotated = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded",
                                                                                             "inFile"] else False
    statusGroup._v_attrs.geneSequences = True if pangenome.status["geneSequences"] in ["Computed", "Loaded",
                                                                                       "inFile"] else False
    statusGroup._v_attrs.genesClustered = True if pangenome.status["genesClustered"] in ["Computed", "Loaded",
                                                                                         "inFile"] else False
    statusGroup._v_attrs.geneFamilySequences = True if pangenome.status["geneFamilySequences"] in ["Computed", "Loaded",
                                                                                                   "inFile"] else False
    statusGroup._v_attrs.NeighborsGraph = True if pangenome.status["neighborsGraph"] in ["Computed", "Loaded",
                                                                                         "inFile"] else False
    statusGroup._v_attrs.Partitionned = True if pangenome.status["partitionned"] in ["Computed", "Loaded",
                                                                                     "inFile"] else False
    statusGroup._v_attrs.defragmented = True if pangenome.status["defragmented"] in ["Computed", "Loaded",
                                                                                     "inFile"] else False
    statusGroup._v_attrs.predictedRGP = True if pangenome.status["predictedRGP"] in ["Computed", "Loaded",
                                                                                     "inFile"] else False
    statusGroup._v_attrs.spots = True if pangenome.status["spots"] in ["Computed", "Loaded", "inFile"] else False
    statusGroup._v_attrs.modules = True if pangenome.status["modules"] in ["Computed", "Loaded", "inFile"] else False

    statusGroup._v_attrs.version = pkg_resources.get_distribution("ppanggolin").version


def writeInfo(pangenome, h5f):
    """ writes information and numbers to be eventually called with the 'info' submodule """
    if "/info" in h5f:
        infoGroup = h5f.root.info
    else:
        infoGroup = h5f.create_group("/", "info", "Informations about the pangenome's content")
    if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]:
        infoGroup._v_attrs.numberOfGenes = len(pangenome.genes)
        infoGroup._v_attrs.numberOfOrganisms = len(pangenome.organisms)
    if pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        infoGroup._v_attrs.numberOfClusters = len(pangenome.geneFamilies)
    if pangenome.status["neighborsGraph"] in ["Computed", "Loaded"]:
        infoGroup._v_attrs.numberOfEdges = len(pangenome.edges)
    if pangenome.status["partitionned"] in ["Computed", "Loaded"]:
        namedPartCounter = Counter()
        subpartCounter = Counter()
        partDistribs = defaultdict(list)
        partSet = set()
        for fam in pangenome.geneFamilies:
            namedPartCounter[fam.namedPartition] += 1
            partDistribs[fam.namedPartition].append(len(fam.organisms) / len(pangenome.organisms))
            if fam.namedPartition == "shell":
                subpartCounter[fam.partition] += 1
            if fam.partition != "S_":
                partSet.add(fam.partition)

        def getmean(arg):
            if len(arg) == 0:
                return 0
            else:
                return round(statistics.mean(arg), 2)

        def getstdev(arg):
            if len(arg) <= 1:
                return 0
            else:
                return round(statistics.stdev(arg), 2)

        def getmax(arg):
            if len(arg) == 0:
                return 0
            else:
                return round(max(arg), 2)

        def getmin(arg):
            if len(arg) == 0:
                return 0
            else:
                return round(min(arg), 2)

        infoGroup._v_attrs.numberOfPersistent = namedPartCounter["persistent"]
        infoGroup._v_attrs.persistentStats = {"min": getmin(partDistribs["persistent"]),
                                              "max": getmax(partDistribs["persistent"]),
                                              "sd": getstdev(partDistribs["persistent"]),
                                              "mean": getmean(partDistribs["persistent"])}
        infoGroup._v_attrs.numberOfShell = namedPartCounter["shell"]
        infoGroup._v_attrs.shellStats = {"min": getmin(partDistribs["shell"]), "max": getmax(partDistribs["shell"]),
                                         "sd": getstdev(partDistribs["shell"]), "mean": getmean(partDistribs["shell"])}
        infoGroup._v_attrs.numberOfCloud = namedPartCounter["cloud"]
        infoGroup._v_attrs.cloudStats = {"min": getmin(partDistribs["cloud"]), "max": getmax(partDistribs["cloud"]),
                                         "sd": getstdev(partDistribs["cloud"]), "mean": getmean(partDistribs["cloud"])}
        infoGroup._v_attrs.numberOfPartitions = len(partSet)
        infoGroup._v_attrs.numberOfSubpartitions = subpartCounter
    if pangenome.status["predictedRGP"] in ["Computed", "Loaded"]:
        infoGroup._v_attrs.numberOfRGP = len(pangenome.regions)
    if pangenome.status["spots"] in ["Computed", "Loaded"]:
        infoGroup._v_attrs.numberOfSpots = len(pangenome.spots)
    if pangenome.status["modules"] in ["Computed", "Loaded"]:
        infoGroup._v_attrs.numberOfModules = len(pangenome.modules)
        infoGroup._v_attrs.numberOfFamiliesInModules = sum([len(mod.families) for mod in pangenome.modules])

    infoGroup._v_attrs.parameters = pangenome.parameters  # saving the pangenome parameters


def updateGeneFamPartition(pangenome, h5f, disable_bar=False):
    logging.getLogger().info("Updating gene families with partition information")
    table = h5f.root.geneFamiliesInfo
    bar = tqdm(range(table.nrows), unit="gene family", disable=disable_bar)
    for row in table:
        row["partition"] = pangenome.getGeneFamily(row["name"].decode()).partition
        row.update()
        bar.update()
    bar.close()


def updateGeneFragments(pangenome, h5f, disable_bar=False):
    """
        updates the annotation table with the fragmentation informations from the defrag pipeline
    """
    logging.getLogger().info("Updating annotations with fragment information")
    table = h5f.root.annotations.genes
    # row = table.row
    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    for row in table:
        if row['gene/type'].decode() == 'CDS':
            row['gene/is_fragment'] = pangenome.getGene(row['gene/ID'].decode()).is_fragment
            row.update()
        bar.update()
    bar.close()
    table.flush()


def ErasePangenome(pangenome, graph=False, geneFamilies=False, partition=False, rgp=False, spots=False, modules=False):
    """ erases tables from a pangenome .h5 file """

    h5f = tables.open_file(pangenome.file, "a")
    statusGroup = h5f.root.status

    if '/edges' in h5f and (graph or geneFamilies):
        logging.getLogger().info("Erasing the formerly computed edges")
        h5f.remove_node("/", "edges")
        statusGroup._v_attrs.NeighborsGraph = False
        pangenome.status["neighborsGraph"] = "No"
    if '/geneFamilies' in h5f and geneFamilies:
        logging.getLogger().info("Erasing the formerly computed gene family to gene associations...")
        h5f.remove_node('/', 'geneFamilies')  # erasing the table, and rewriting a new one.
        pangenome.status["defragmented"] = "No"
        pangenome.status["genesClustered"] = "No"
        statusGroup._v_attrs.defragmented = False
        statusGroup._v_attrs.genesClustered = False
    if '/geneFamiliesInfo' in h5f and geneFamilies:
        logging.getLogger().info("Erasing the formerly computed gene family representative sequences...")
        h5f.remove_node('/', 'geneFamiliesInfo')  # erasing the table, and rewriting a new one.
        pangenome.status["partitionned"] = "No"
        pangenome.status["geneFamilySequences"] = "No"
        statusGroup._v_attrs.geneFamilySequences = False
        statusGroup._v_attrs.Partitionned = False
    if '/RGP' in h5f and (geneFamilies or partition or rgp):
        logging.getLogger().info("Erasing the formerly computer RGP...")
        pangenome.status["predictedRGP"] = "No"
        statusGroup._v_attrs.predictedRGP = False
        h5f.remove_node("/", "RGP")
    if '/spots' in h5f and (geneFamilies or partition or rgp or spots):
        logging.getLogger().info("Erasing the formerly computed spots...")
        pangenome.status["spots"] = "No"
        statusGroup._v_attrs.spots = False
        h5f.remove_node("/", "spots")

    if '/modules' in h5f and (geneFamilies or partition or modules):
        logging.getLogger().info("Erasing the formerly computed modules...")
        pangenome.status["modules"] = "No"
        statusGroup._v_attrs.modules = False
        h5f.remove_node("/", "modules")

    h5f.close()


def writePangenome(pangenome, filename, force, disable_bar=False):
    """
        Writes or updates a pangenome file
        pangenome is the corresponding pangenome object, filename the h5 file and status what has been modified.
    """

    if pangenome.status["genomesAnnotated"] == "Computed":
        compressionFilter = tables.Filters(complevel=1, shuffle=True, bitshuffle=True, complib='blosc:zstd')
        h5f = tables.open_file(filename, "w", filters=compressionFilter)
        logging.getLogger().info("Writing genome annotations...")
        writeAnnotations(pangenome, h5f)
        pangenome.status["genomesAnnotated"] = "Loaded"
        h5f.close()
    elif pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"]:
        pass
    else:  # if the pangenome is not Computed not Loaded, it's probably not really in a good state ( or something new was coded).
        raise NotImplementedError(
            "Something REALLY unexpected and unplanned for happened here. Please post an issue on github with what you did to reach this error.")

    # from there, appending to existing file.
    h5f = tables.open_file(filename, "a")

    if pangenome.status["geneSequences"] == "Computed":
        logging.getLogger().info("writing the protein coding gene dna sequences")
        writeGeneSequences(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["geneSequences"] = "Loaded"

    if pangenome.status["genesClustered"] == "Computed":
        logging.getLogger().info("Writing gene families and gene associations...")
        writeGeneFamilies(pangenome, h5f, force, disable_bar=disable_bar)
        logging.getLogger().info("Writing gene families information...")
        writeGeneFamInfo(pangenome, h5f, force, disable_bar=disable_bar)
        if pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"] and pangenome.status[
            "defragmented"] == "Computed":
            # if the annotations have not been computed in this run, and there has been a clustering with defragmentation, then the annotations can be updated
            updateGeneFragments(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["genesClustered"] = "Loaded"
    if pangenome.status["neighborsGraph"] == "Computed":
        logging.getLogger().info("Writing the edges...")
        writeGraph(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status["neighborsGraph"] = "Loaded"

    if pangenome.status["partitionned"] == "Computed" and pangenome.status["genesClustered"] in ["Loaded",
                                                                                                 "inFile"]:  # otherwise it's been written already.
        updateGeneFamPartition(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["partitionned"] = "Loaded"

    if pangenome.status['predictedRGP'] == "Computed":
        logging.getLogger().info("Writing Regions of Genomic Plasticity...")
        writeRGP(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status['predictedRGP'] = "Loaded"

    if pangenome.status["spots"] == "Computed":
        logging.getLogger().info("Writing Spots of Insertion...")
        writeSpots(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status['spots'] = "Loaded"

    if pangenome.status["modules"] == "Computed":
        logging.getLogger().info("Writing Modules...")
        writeModules(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status["modules"] = "Loaded"

    writeStatus(pangenome, h5f)
    writeInfo(pangenome, h5f)

    h5f.close()
    logging.getLogger().info(f"Done writing the pangenome. It is in file : {filename}")
