#!/usr/bin/env python3

# default libraries
import logging
from collections import Counter, defaultdict
import statistics
from typing import Tuple, Union
from importlib.metadata import distribution
import os

# installed libraries
from tqdm import tqdm
import tables
from gmpy2 import popcount

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats.writeAnnotations import write_annotations, write_gene_sequences
from ppanggolin.formats.writeMetadata import (
    write_metadata,
    erase_metadata,
    write_metadata_status,
)
from ppanggolin.genome import Feature, Gene
from ppanggolin.formats.readBinaries import read_genedata, Genedata


def getmean(arg: iter) -> float:
    """Compute the mean of arguments if exist 0 else

    :param arg: list of values

    :return: return the mean
    """
    return 0 if len(arg) == 0 else round(statistics.mean(arg), 2)


def getstdev(arg: iter) -> float:
    """Compute the standard deviation of arguments if exist 0 else

    :param arg: list of values

    :return: return the sd
    """
    return 0 if len(arg) <= 1 else round(statistics.stdev(arg), 2)


def getmax(arg: iter) -> float:
    """Get the maximum of arguments if exist 0 else

    :param arg: list of values

    :return: return the maximum
    """
    return 0 if len(arg) == 0 else round(max(arg), 2)


def getmin(arg: iter) -> float:
    """Get the minimum of arguments if exist 0 else

    :param arg: list of values

    :return: return the minimum
    """
    return 0 if len(arg) == 0 else round(min(arg), 2)


def gene_fam_desc(
    max_name_len: int, max_sequence_length: int, max_part_len: int
) -> dict:
    """
    Create a formatted table for gene families description

    :param max_name_len: Maximum size of gene family name
    :param max_sequence_length: Maximum size of gene family representing gene sequences
    :param max_part_len: Maximum size of gene family partition

    :return: Formatted table
    """
    return {
        "name": tables.StringCol(itemsize=max_name_len),
        "protein": tables.StringCol(itemsize=max_sequence_length),
        "partition": tables.StringCol(itemsize=max_part_len),
    }


def get_gene_fam_len(pangenome: Pangenome) -> Tuple[int, int, int]:
    """
    Get maximum size of gene families information

    :param pangenome: Pangenome with gene families computed

    :return: Maximum size of each element
    """
    max_gene_fam_name_len = 1
    max_gene_fam_seq_len = 1
    max_part_len = 3
    for genefam in pangenome.gene_families:
        if len(genefam.sequence) > max_gene_fam_seq_len:
            max_gene_fam_seq_len = len(genefam.sequence)
        if len(genefam.name) > max_gene_fam_name_len:
            max_gene_fam_name_len = len(genefam.name)
        if len(genefam.partition) > max_part_len:
            max_part_len = len(genefam.partition)
    return max_gene_fam_name_len, max_gene_fam_seq_len, max_part_len


def write_gene_fam_info(
    pangenome: Pangenome,
    h5f: tables.File,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Writing a table containing the protein sequences of each family

    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param force: force to write information if precedent information exist
    :param disable_bar: Disable progress bar
    """
    if "/geneFamiliesInfo" in h5f and force is True:
        logging.getLogger("PPanGGOLiN").info(
            "Erasing the formerly computed gene family representative sequences..."
        )
        h5f.remove_node(
            "/", "geneFamiliesInfo"
        )  # erasing the table, and rewriting a new one.
    gene_fam_seq = h5f.create_table(
        "/",
        "geneFamiliesInfo",
        gene_fam_desc(*get_gene_fam_len(pangenome)),
        expectedrows=pangenome.number_of_gene_families,
    )

    row = gene_fam_seq.row
    for fam in tqdm(
        pangenome.gene_families,
        total=pangenome.number_of_gene_families,
        unit="gene family",
        disable=disable_bar,
    ):
        row["name"] = fam.name
        row["protein"] = fam.sequence
        row["partition"] = fam.partition
        row.append()
    gene_fam_seq.flush()


def gene_to_fam_desc(gene_fam_name_len: int, gene_id_len: int) -> dict:
    """
    Create a formatted table for gene in gene families information

    :param gene_fam_name_len: Maximum size of gene family names
    :param gene_id_len: Maximum size of gene identifier

    :return: formatted table
    """
    return {
        "geneFam": tables.StringCol(itemsize=gene_fam_name_len),
        "gene": tables.StringCol(itemsize=gene_id_len),
    }


def get_gene_to_fam_len(pangenome: Pangenome):
    """
    Get maximum size of gene in gene families information

    :param pangenome: Pangenome with gene families computed

    :return: Maximum size of each element
    """
    max_gene_fam_name = 1
    max_gene_id = 1
    for family in pangenome.gene_families:
        if len(family.name) > max_gene_fam_name:
            max_gene_fam_name = len(family.name)
        for gene in family.genes:
            if len(gene.ID) > max_gene_id:
                max_gene_id = len(gene.ID)
    return max_gene_fam_name, max_gene_id


def write_gene_families(
    pangenome: Pangenome,
    h5f: tables.File,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Function writing all the pangenome gene families

    :param pangenome: pangenome with gene families computed
    :param h5f: HDF5 file to save pangenome with gene families
    :param force: Force to write gene families in hdf5 file if there is already gene families
    :param disable_bar: Disable progress bar
    """
    if "/geneFamilies" in h5f and force is True:
        logging.getLogger("PPanGGOLiN").info(
            "Erasing the formerly computed gene family to gene associations..."
        )
        h5f.remove_node(
            "/", "geneFamilies"
        )  # erasing the table, and rewriting a new one.
    gene_families = h5f.create_table(
        "/", "geneFamilies", gene_to_fam_desc(*get_gene_to_fam_len(pangenome))
    )
    gene_row = gene_families.row
    for family in tqdm(
        pangenome.gene_families,
        total=pangenome.number_of_gene_families,
        unit="gene family",
        disable=disable_bar,
    ):
        for gene in family.genes:
            gene_row["gene"] = gene.ID
            gene_row["geneFam"] = family.name
            gene_row.append()
    gene_families.flush()


def graph_desc(max_gene_id_len):
    """
    Create a formatted table for pangenome graph

    :param max_gene_id_len: Maximum size of gene id

    :return: formatted table
    """
    return {
        "geneTarget": tables.StringCol(itemsize=max_gene_id_len),
        "geneSource": tables.StringCol(itemsize=max_gene_id_len),
    }


def get_gene_id_len(pangenome: Pangenome) -> int:
    """
    Get maximum size of gene id in pangenome graph

    :param pangenome: Pangenome with graph computed

    :return: Maximum size of gene id
    """
    max_gene_len = 1
    for gene in pangenome.genes:
        if len(gene.ID) > max_gene_len:
            max_gene_len = len(gene.ID)
    return max_gene_len


def write_graph(
    pangenome: Pangenome,
    h5f: tables.File,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Function writing the pangenome graph

    :param pangenome: pangenome with graph computed
    :param h5f: HDF5 file to save pangenome graph
    :param force: Force to write graph in hdf5 file if there is already one
    :param disable_bar: Disable progress bar
    """
    # TODO if we want to be able to read the graph without reading the annotations (because it's one of the most time
    #  consumming parts to read), it might be good to add the organism name in the table here.
    #  for now, forcing the read of annotations.
    if "/edges" in h5f and force is True:
        logging.getLogger("PPanGGOLiN").info("Erasing the formerly computed edges")
        h5f.remove_node("/", "edges")
    edge_table = h5f.create_table(
        "/",
        "edges",
        graph_desc(get_gene_id_len(pangenome)),
        expectedrows=pangenome.number_of_edges,
    )
    edge_row = edge_table.row
    for edge in tqdm(
        pangenome.edges,
        total=pangenome.number_of_edges,
        unit="edge",
        disable=disable_bar,
    ):
        for gene1, gene2 in edge.gene_pairs:
            edge_row["geneTarget"] = gene1.ID
            edge_row["geneSource"] = gene2.ID
            edge_row.append()
    edge_table.flush()


def rgp_desc(max_rgp_len, max_gene_len):
    """
    Create a formatted table for region of genomic plasticity

    :param max_rgp_len: Maximum size of RGP
    :param max_gene_len: Maximum sizez of gene

    :return: formatted table
    """
    return {
        "RGP": tables.StringCol(itemsize=max_rgp_len),
        "gene": tables.StringCol(itemsize=max_gene_len),
        "score": tables.UInt32Col(),
    }


def get_rgp_len(pangenome: Pangenome) -> Tuple[int, int]:
    """
    Get maximum size of region of genomic plasticity and gene

    :param pangenome: Pangenome with gene families computed

    :return: Maximum size of each element
    """
    max_gene_len = 1
    max_rgp_len = 1
    for region in pangenome.regions:
        for gene in region.genes:
            if len(gene.ID) > max_gene_len:
                max_gene_len = len(gene.ID)
        if len(region.name) > max_rgp_len:
            max_rgp_len = len(region.name)
    return max_rgp_len, max_gene_len


def write_rgp(
    pangenome: Pangenome,
    h5f: tables.File,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Function writing all the region of genomic plasticity in pangenome

    :param pangenome: pangenome with RGP computed
    :param h5f: HDF5 file to save pangenome with RGP
    :param force: Force to write gene families in hdf5 file if there is already RGP
    :param disable_bar: Disable progress bar
    """
    if "/RGP" in h5f and force is True:
        logging.getLogger("PPanGGOLiN").info("Erasing the formerly computer RGP")
        h5f.remove_node("/", "RGP")

    rgp_table = h5f.create_table(
        "/",
        "RGP",
        rgp_desc(*get_rgp_len(pangenome)),
        expectedrows=sum([len(region) for region in pangenome.regions]),
    )
    rgp_row = rgp_table.row
    for region in tqdm(
        pangenome.regions,
        total=pangenome.number_of_rgp,
        unit="region",
        disable=disable_bar,
    ):
        for gene in region.genes:
            rgp_row["RGP"] = region.name
            rgp_row["gene"] = gene.ID
            rgp_row["score"] = region.score
            rgp_row.append()
    rgp_table.flush()


def spot_desc(max_rgp_len):
    """
    Create a formatted table for hotspot

    :param max_rgp_len: Maximum size of RGP

    :return: formatted table
    """
    return {"spot": tables.UInt32Col(), "RGP": tables.StringCol(itemsize=max_rgp_len)}


def get_spot_desc(pangenome: Pangenome) -> int:
    """
    Get maximum size of region of genomic plasticity in hotspot

    :param pangenome: Pangenome with gene families computed

    :return: Maximum size of each element
    """
    max_rgp_len = 1
    for spot in pangenome.spots:
        for region in spot.regions:
            if len(region.name) > max_rgp_len:
                max_rgp_len = len(region.name)
    return max_rgp_len


def write_spots(
    pangenome: Pangenome,
    h5f: tables.File,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Function writing all the pangenome hotspot

    :param pangenome: pangenome with spot computed
    :param h5f: HDF5 file to save pangenome with spot
    :param force: Force to write gene families in hdf5 file if there is already spot
    :param disable_bar: Disable progress bar
    """
    if "/spots" in h5f and force is True:
        logging.getLogger("PPanGGOLiN").info("Erasing the formerly computed spots")
        h5f.remove_node("/", "spots")

    spot_table = h5f.create_table(
        "/",
        "spots",
        spot_desc(get_spot_desc(pangenome)),
        expectedrows=sum([len(spot) for spot in pangenome.spots]),
    )
    spot_row = spot_table.row
    for spot in tqdm(
        pangenome.spots,
        total=pangenome.number_of_spots,
        unit="spot",
        disable=disable_bar,
    ):
        for region in spot.regions:
            spot_row["spot"] = spot.ID
            spot_row["RGP"] = region.name
            spot_row.append()
    spot_table.flush()


def mod_desc(gene_fam_name_len):
    """
    Create a formatted table for hotspot

    :param gene_fam_name_len: Maximum size of gene families name

    :return: formatted table
    """
    return {
        "geneFam": tables.StringCol(itemsize=gene_fam_name_len),
        "module": tables.UInt32Col(),
    }


def get_mod_desc(pangenome: Pangenome) -> int:
    """
    Get maximum size of gene families name in modules

    :param pangenome: Pangenome with modules computed

    :return: Maximum size of each element
    """
    max_fam_len = 1
    for mod in pangenome.modules:
        for fam in mod.families:
            if len(fam.name) > max_fam_len:
                max_fam_len = len(fam.name)
    return max_fam_len


def write_modules(
    pangenome: Pangenome,
    h5f: tables.File,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Function writing all the pangenome modules

    :param pangenome: pangenome with spot computed
    :param h5f: HDF5 file to save pangenome with spot
    :param force: Force to write gene families in hdf5 file if there is already spot
    :param disable_bar: Disable progress bar
    """
    if "/modules" in h5f and force is True:
        logging.getLogger("PPanGGOLiN").info("Erasing the formerly computed modules")
        h5f.remove_node("/", "modules")

    mod_table = h5f.create_table(
        "/",
        "modules",
        mod_desc(get_mod_desc(pangenome)),
        expectedrows=sum([len(mod) for mod in pangenome.modules]),
    )
    mod_row = mod_table.row

    for mod in tqdm(
        pangenome.modules,
        total=pangenome.number_of_modules,
        unit="modules",
        disable=disable_bar,
    ):
        for fam in mod.families:
            mod_row["geneFam"] = fam.name
            mod_row["module"] = mod.ID
            mod_row.append()
    mod_table.flush()

    write_info_modules(pangenome, h5f)


def write_status(pangenome: Pangenome, h5f: tables.File):
    """
    Write pangenome status in HDF5 file

    :param pangenome: Pangenome object
    :param h5f: Pangenome file
    """
    if "/status" in h5f:  # if statuses are already written
        status_group = h5f.root.status
    else:  # else create the status group.
        status_group = h5f.create_group(
            "/", "status", "Statuses of the pangenome content"
        )
    status_group._v_attrs.genomesAnnotated = (
        True
        if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.geneSequences = (
        True
        if pangenome.status["geneSequences"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.genesClustered = (
        True
        if pangenome.status["genesClustered"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.geneFamilySequences = (
        True
        if pangenome.status["geneFamilySequences"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.NeighborsGraph = (
        True
        if pangenome.status["neighborsGraph"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.Partitioned = (
        True
        if pangenome.status["partitioned"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.defragmented = (
        True
        if pangenome.status["defragmented"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.predictedRGP = (
        True
        if pangenome.status["predictedRGP"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.spots = (
        True if pangenome.status["spots"] in ["Computed", "Loaded", "inFile"] else False
    )
    status_group._v_attrs.modules = (
        True
        if pangenome.status["modules"] in ["Computed", "Loaded", "inFile"]
        else False
    )
    status_group._v_attrs.metadata = write_metadata_status(pangenome, h5f, status_group)
    status_group._v_attrs.version = distribution("ppanggolin").version


def write_info(pangenome: Pangenome, h5f: tables.File):
    """
    Writes information and numbers to be eventually called with the 'info' submodule

    :param pangenome: Pangenome object with some information computed
    :param h5f: Pangenome file to save information
    """

    if "/info" in h5f:
        info_group = h5f.root.info
    else:
        info_group = h5f.create_group(
            "/", "info", "Information about the pangenome content"
        )
    if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfGenes = pangenome.number_of_genes
        info_group._v_attrs.numberOfGenomes = pangenome.number_of_organisms
    if pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfClusters = pangenome.number_of_gene_families
    if pangenome.status["neighborsGraph"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfEdges = pangenome.number_of_edges
    if pangenome.status["partitioned"] in ["Computed", "Loaded"]:
        named_part_counter = Counter()
        subpart_counter = Counter()
        part_distribs = defaultdict(list)
        part_set = set()
        for fam in pangenome.gene_families:
            named_part_counter[fam.named_partition] += 1
            part_distribs[fam.named_partition].append(
                fam.number_of_organisms / pangenome.number_of_organisms
            )
            if fam.named_partition == "shell":
                subpart_counter[fam.partition] += 1
            if fam.partition != "S_":
                part_set.add(fam.partition)

        info_group._v_attrs.numberOfPersistent = named_part_counter["persistent"]
        info_group._v_attrs.persistentStats = {
            "min_genomes_frequency": getmin(part_distribs["persistent"]),
            "max_genomes_frequency": getmax(part_distribs["persistent"]),
            "sd_genomes_frequency": getstdev(part_distribs["persistent"]),
            "mean_genomes_frequency": getmean(part_distribs["persistent"]),
        }

        info_group._v_attrs.numberOfShell = named_part_counter["shell"]
        info_group._v_attrs.shellStats = {
            "min_genomes_frequency": getmin(part_distribs["shell"]),
            "max_genomes_frequency": getmax(part_distribs["shell"]),
            "sd_genomes_frequency": getstdev(part_distribs["shell"]),
            "mean_genomes_frequency": getmean(part_distribs["shell"]),
        }

        info_group._v_attrs.numberOfCloud = named_part_counter["cloud"]
        info_group._v_attrs.cloudStats = {
            "min_genomes_frequency": getmin(part_distribs["cloud"]),
            "max_genomes_frequency": getmax(part_distribs["cloud"]),
            "sd_genomes_frequency": getstdev(part_distribs["cloud"]),
            "mean_genomes_frequency": getmean(part_distribs["cloud"]),
        }

        info_group._v_attrs.numberOfPartitions = len(part_set)
        info_group._v_attrs.numberOfSubpartitions = subpart_counter

    if pangenome.status["predictedRGP"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfRGP = pangenome.number_of_rgp

    if pangenome.status["spots"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfSpots = pangenome.number_of_spots

    if pangenome.status["modules"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfModules = pangenome.number_of_modules
        info_group._v_attrs.numberOfFamiliesInModules = sum(
            [len(mod) for mod in pangenome.modules]
        )

    info_group._v_attrs.parameters = (
        pangenome.parameters
    )  # saving the pangenome parameters


def write_info_modules(pangenome: Pangenome, h5f: tables.File):
    """
    Writes information about modules

    :param pangenome: Pangenome object with some information computed
    :param h5f: Pangenome file to save information
    """

    def part_spec(part: str) -> list:
        """
        Get the list of module for a specific partition of pangenome

        :param part: pangenome partition name

        :return: list of module specific to partition
        """
        pangenome.compute_mod_bitarrays(part)
        return [popcount(module.bitarray) for module in pangenome.modules]

    if "/info" not in h5f:
        write_info(pangenome, h5f)
    info_group = h5f.root.info

    mod_fam = [len(module) for module in pangenome.modules]
    sum_mod_fam = sum(mod_fam)

    info_group._v_attrs.StatOfFamiliesInModules = {
        "min": getmin(mod_fam),
        "max": getmax(mod_fam),
        "sd": getstdev(mod_fam),
        "mean": getmean(mod_fam),
    }

    spec_pers = part_spec(part="persistent")
    spec_shell = part_spec(part="shell")
    spec_cloud = part_spec(part="cloud")

    info_group._v_attrs.PersistentSpecInModules = {
        "percent": (
            round((sum(spec_pers) / sum_mod_fam) * 100, 2) if sum_mod_fam > 0 else 0
        ),
        "min": getmin(spec_pers),
        "max": getmax(spec_pers),
        "sd": getstdev(spec_pers),
        "mean": getmean(spec_pers),
    }

    info_group._v_attrs.ShellSpecInModules = {
        "percent": (
            round((sum(spec_shell) / sum_mod_fam) * 100, 2) if sum_mod_fam > 0 else 0
        ),
        "min": getmin(spec_shell),
        "max": getmax(spec_shell),
        "sd": getstdev(spec_shell),
        "mean": getmean(spec_shell),
    }

    info_group._v_attrs.CloudSpecInModules = {
        "percent": (
            round((sum(spec_cloud) / sum_mod_fam) * 100, 2) if sum_mod_fam > 0 else 0
        ),
        "min": getmin(spec_cloud),
        "max": getmax(spec_cloud),
        "sd": getstdev(spec_cloud),
        "mean": getmean(spec_cloud),
    }


def update_gene_fam_partition(
    pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False
):
    """
    Update the gene families table with partition information

    :param pangenome: Partitioned pangenome
    :param h5f: HDF5 file with gene families
    :param disable_bar: Allow to disable progress bar
    """
    logging.getLogger("PPanGGOLiN").info(
        "Updating gene families with partition information"
    )
    table = h5f.root.geneFamiliesInfo
    for row in tqdm(table, total=table.nrows, unit="gene family", disable=disable_bar):
        row["partition"] = pangenome.get_gene_family(row["name"].decode()).partition
        row.update()


def update_gene_fragments(
    pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False
):
    """
    Updates the annotation table with the fragmentation information from the defrag pipeline

    :param pangenome: Annotated pangenome
    :param h5f: HDF5 pangenome file
    :param disable_bar: Allow to disable progress bar
    """
    logging.getLogger("PPanGGOLiN").info(
        "Updating annotations with fragment information"
    )
    genedataid2genedata = read_genedata(h5f)

    table = h5f.root.annotations.genes
    for row in tqdm(table, total=table.nrows, unit="gene", disable=disable_bar):
        genedata_id = row["genedata_id"]
        if genedataid2genedata[genedata_id].gene_type == "CDS":
            row["is_fragment"] = pangenome.get_gene(row["ID"].decode()).is_fragment
            row.update()
    table.flush()


def erase_pangenome(
    pangenome: Pangenome,
    graph: bool = False,
    gene_families: bool = False,
    partition: bool = False,
    rgp: bool = False,
    spots: bool = False,
    modules: bool = False,
    metadata: bool = False,
    metatype: str = None,
    source: str = None,
):
    """
    Erases tables from a pangenome .h5 file

    :param pangenome: Pangenome
    :param graph: remove graph information
    :param gene_families: remove gene families information
    :param partition: remove partition information
    :param rgp: remove rgp information
    :param spots: remove spots information
    :param modules: remove modules information
    :param metadata: remove metadata information
    :param metatype:
    :param source:
    """

    if metadata and (metatype is None or source is None):
        raise ValueError("To erase metadata. You should provide metatype and source")

    with tables.open_file(pangenome.file, "a") as h5f:
        status_group = h5f.root.status
        info_group = h5f.root.info

        if "/edges" in h5f and (graph or gene_families):
            logging.getLogger("PPanGGOLiN").info("Erasing the formerly computed edges")
            h5f.remove_node("/", "edges")
            status_group._v_attrs.NeighborsGraph = False
            pangenome.status["neighborsGraph"] = "No"
            h5f.del_node_attr(info_group, "numberOfEdges")
        if "/geneFamilies" in h5f and gene_families:
            logging.getLogger("PPanGGOLiN").info(
                "Erasing the formerly computed gene family to gene associations..."
            )
            h5f.remove_node(
                "/", "geneFamilies"
            )  # erasing the table, and rewriting a new one.
            pangenome.status["defragmented"] = "No"
            pangenome.status["genesClustered"] = "No"
            status_group._v_attrs.defragmented = False
            status_group._v_attrs.genesClustered = False

            h5f.del_node_attr(info_group, "numberOfClusters")

        if "/geneFamiliesInfo" in h5f and gene_families:
            logging.getLogger("PPanGGOLiN").info(
                "Erasing the formerly computed gene family representative sequences..."
            )
            h5f.remove_node(
                "/", "geneFamiliesInfo"
            )  # erasing the table, and rewriting a new one.
            pangenome.status["geneFamilySequences"] = "No"
            status_group._v_attrs.geneFamilySequences = False
            if partition:
                logging.getLogger("PPanGGOLiN").info("Erasing former partitions...")
                pangenome.status["partitioned"] = "No"
                status_group._v_attrs.Partitioned = False
                if "Partitioned" in status_group._v_attrs._f_list():
                    status_group._v_attrs.Partitioned = False

                h5f.del_node_attr(info_group, "numberOfPersistent")
                h5f.del_node_attr(info_group, "persistentStats")
                h5f.del_node_attr(info_group, "numberOfShell")
                h5f.del_node_attr(info_group, "shellStats")
                h5f.del_node_attr(info_group, "numberOfCloud")
                h5f.del_node_attr(info_group, "cloudStats")
                h5f.del_node_attr(info_group, "numberOfPartitions")
                h5f.del_node_attr(info_group, "numberOfSubpartitions")

        if "/RGP" in h5f and (gene_families or partition or rgp):
            logging.getLogger("PPanGGOLiN").info("Erasing the formerly computer RGP...")
            pangenome.status["predictedRGP"] = "No"
            status_group._v_attrs.predictedRGP = False
            h5f.remove_node("/", "RGP")

            h5f.del_node_attr(info_group, "numberOfRGP")

        if "/spots" in h5f and (gene_families or partition or rgp or spots):
            logging.getLogger("PPanGGOLiN").info(
                "Erasing the formerly computed spots..."
            )
            pangenome.status["spots"] = "No"
            status_group._v_attrs.spots = False
            h5f.remove_node("/", "spots")

            h5f.del_node_attr(info_group, "numberOfSpots")

        if "/modules" in h5f and (gene_families or partition or modules):
            logging.getLogger("PPanGGOLiN").info(
                "Erasing the formerly computed modules..."
            )
            pangenome.status["modules"] = "No"
            status_group._v_attrs.modules = False
            h5f.remove_node("/", "modules")

            h5f.del_node_attr(info_group, "numberOfModules")
            h5f.del_node_attr(info_group, "numberOfFamiliesInModules")
            for info in [
                "CloudSpecInModules",
                "PersistentSpecInModules",
                "ShellSpecInModules",
                "numberOfFamiliesInModules",
                "StatOfFamiliesInModules",
            ]:
                if info in info_group._v_attrs._f_list():
                    h5f.del_node_attr(info_group, info)

        if "/metadata/" in h5f and metadata:
            erase_metadata(pangenome, h5f, status_group, metatype, source)


def write_pangenome(
    pangenome: Pangenome, filename, force: bool = False, disable_bar: bool = False
):
    """
    Writes or updates a pangenome file

    :param pangenome: pangenome object
    :param filename: HDF5 file to save pangenome
    :param force: force to write on pangenome if information already exist
    :param disable_bar: Allow to disable progress bar
    """
    try:
        assert pangenome.status["genomesAnnotated"] in ["Computed", "Loaded", "inFile"]
    except AssertionError:
        raise AssertionError(
            "Something REALLY unexpected and unplanned for happened here. "
            "Please post an issue on github with what you did to reach this error."
        )

    if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded", "inFile"]:
        if pangenome.status["genomesAnnotated"] == "Computed":
            complevel_default = 6
            try:
                complevel = int(
                    os.getenv("PPANGGOLIN_HDF5_COMPRESSION_LEVEL", complevel_default)
                )
            except ValueError:
                logging.getLogger("PPanGGOLiN").info(
                    f"Invalid PPANGGOLIN_HDF5_COMPRESSION_LEVEL, using default value {complevel_default}"
                )
                complevel = complevel_default

            logging.getLogger("PPanGGOLiN").debug(
                f"Using compression level {complevel} for the pangenome HDF5 file."
            )
            compression_filter = tables.Filters(
                complevel=complevel,
                shuffle=False,
                bitshuffle=False,
                complib="blosc2:zstd",
            )

            h5f = tables.open_file(filename, "w", filters=compression_filter)
            logging.getLogger("PPanGGOLiN").info("Writing genome annotations...")

            write_annotations(pangenome, h5f, disable_bar=disable_bar)

            pangenome.status["genomesAnnotated"] = "Loaded"
            h5f.close()

    # from there, appending to existing file
    h5f = tables.open_file(filename, "a")

    if pangenome.status["geneSequences"] == "Computed":
        logging.getLogger("PPanGGOLiN").info(
            "writing the protein coding gene dna sequences in pangenome..."
        )
        write_gene_sequences(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["geneSequences"] = "Loaded"

    if pangenome.status["genesClustered"] == "Computed":
        logging.getLogger("PPanGGOLiN").info(
            "Writing gene families and gene associations in pangenome..."
        )
        write_gene_families(pangenome, h5f, force, disable_bar=disable_bar)
        logging.getLogger("PPanGGOLiN").info(
            "Writing gene families information in pangenome..."
        )
        write_gene_fam_info(pangenome, h5f, force, disable_bar=disable_bar)
        if (
            pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"]
            and pangenome.status["defragmented"] == "Computed"
        ):
            # if the annotations have not been computed in this run,
            # and there has been a clustering with defragmentation, then the annotations can be updated
            update_gene_fragments(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["genesClustered"] = "Loaded"
    if pangenome.status["neighborsGraph"] == "Computed":
        logging.getLogger("PPanGGOLiN").info(
            "Writing the edges of neighbors graph in pangenome..."
        )
        write_graph(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status["neighborsGraph"] = "Loaded"
    if pangenome.status["partitioned"] == "Computed" and pangenome.status[
        "genesClustered"
    ] in [
        "Loaded",
        "inFile",
    ]:  # otherwise, it's been written already.
        update_gene_fam_partition(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["partitioned"] = "Loaded"

    if pangenome.status["predictedRGP"] == "Computed":
        logging.getLogger("PPanGGOLiN").info(
            "Writing Regions of Genomic Plasticity in pangenome..."
        )
        write_rgp(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status["predictedRGP"] = "Loaded"

    if pangenome.status["spots"] == "Computed":
        logging.getLogger("PPanGGOLiN").info(
            "Writing Spots of Insertion in pangenome..."
        )
        write_spots(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status["spots"] = "Loaded"

    if pangenome.status["modules"] == "Computed":
        logging.getLogger("PPanGGOLiN").info("Writing Modules in pangenome...")
        write_modules(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status["modules"] = "Loaded"

    write_metadata(pangenome, h5f, disable_bar)

    write_status(pangenome, h5f)
    write_info(pangenome, h5f)

    h5f.close()
    logging.getLogger("PPanGGOLiN").info(
        f"Done writing the pangenome. It is in file : {filename}"
    )
