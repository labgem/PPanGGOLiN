#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import sys

# installed libraries
from typing import TextIO

from tables import Table
from tqdm import tqdm
import tables

# local libraries
from ppanggolin.genome import Organism, Gene, RNA
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Spot, Module


def get_number_of_organisms(pangenome: Pangenome) -> int:
    """ Standalone function to get the number of organisms in a pangenome

    :param pangenome: Annotated pangenome

    :return: Number of organisms in the pangenome
    """
    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")
    annotations = h5f.root.annotations

    table = annotations.genes
    org_set = set()
    for org in read_chunks(table, column="organism"):
        org_set.add(org)
    h5f.close()
    return len(org_set)


def get_status(pangenome: Pangenome, pangenome_file: str):
    """
    Checks which elements are already present in the file.

    :param pangenome: Blank pangenome
    :param pangenome_file: path to the pangenome file
    """
    h5f = tables.open_file(pangenome_file, "r")
    logging.getLogger().info("Getting the current pangenome status")
    status_group = h5f.root.status
    if status_group._v_attrs.genomesAnnotated:
        pangenome.status["genomesAnnotated"] = "inFile"
    if status_group._v_attrs.genesClustered:
        pangenome.status["genesClustered"] = "inFile"
    if status_group._v_attrs.geneSequences:
        pangenome.status["geneSequences"] = "inFile"
    if status_group._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "inFile"
    if status_group._v_attrs.NeighborsGraph:
        pangenome.status["neighborsGraph"] = "inFile"

    if 'Partitionned' in status_group._v_attrs._f_list():
        # Partitionned keep working with older version
        h5f.close()
        h5f = tables.open_file(pangenome_file, "a")
        status_group = h5f.root.status
        if status_group._v_attrs.Partitionned:
            status_group._v_attrs.Partitioned = True
        else:
            status_group._v_attrs.Partitioned = False
        del status_group._v_attrs.Partitionned

    if status_group._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "inFile"

    if hasattr(status_group._v_attrs, "predictedRGP") and status_group._v_attrs.predictedRGP:
        pangenome.status["predictedRGP"] = "inFile"

    if hasattr(status_group._v_attrs, "spots") and status_group._v_attrs.spots:
        pangenome.status["spots"] = "inFile"

    if hasattr(status_group._v_attrs, "modules") and status_group._v_attrs.modules:
        pangenome.status["modules"] = "inFile"

    if "/info" in h5f:
        info_group = h5f.root.info
        pangenome.parameters = info_group._v_attrs.parameters
    h5f.close()


def read_chunks(table: Table, column: str = None, chunk: int = 10000):
    """
    Reading entirely the provided table (or column if specified) chunk per chunk to limit RAM usage.

    :param table:
    :param column:
    :param chunk:
    """
    for i in range(0, table.nrows, chunk):
        for row in table.read(start=i, stop=i + chunk, field=column):
            yield row


def read_sequences(h5f: tables.File) -> dict:
    """
    Reads the sequences table and returns a seqid2seq dictionnary
    :param h5f: the hdf5 file handler
    :return: dictionnary linking sequences to the seq identifier
    """
    table = h5f.root.sequences
    seqid2seq = {}
    for row in read_chunks(table,chunk=20000):
        seqid2seq[row["seqid"].decode()] = row['dna'].decode()
    return seqid2seq

def get_gene_sequences_from_file(filename: str, file_obj: TextIO, list_cds: iter = None, add: str = '',
                                 disable_bar: bool = False):
    """
    Writes the CDS sequences of the Pangenome object to a File object that can be filtered or not by a list of CDS,
    and adds the eventual str 'add' in front of the identifiers. Loads the sequences from a .h5 pangenome file.

    :param filename: Name of the pangenome file
    :param file_obj: Name of the output file
    :param list_cds: An iterable object of CDS
    :param add: Add a prefix to sequence header
    :param disable_bar: Prevent to print disable progress bar
    """
    logging.getLogger().info(f"Extracting and writing CDS sequences from a {filename} file to a fasta file...")
    h5f = tables.open_file(filename, "r", driver_core_backing_store=0)
    table = h5f.root.geneSequences
    list_cds = set(list_cds) if list_cds is not None else None
    seqid2seq = read_sequences(h5f)
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene", disable=disable_bar):
        # Read the table chunk per chunk otherwise RAM dies on big pangenomes
        name_cds = row["gene"].decode()
        if row["type"] == b"CDS" and (list_cds is None or name_cds in list_cds):
            file_obj.write('>' + add + name_cds + "\n")
            file_obj.write(seqid2seq[row["seqid"].decode()] + "\n")
    file_obj.flush()
    h5f.close()


def launch_read_organism(args) -> None:
    """
    Allow to launch read organism in multiprocessing

    :param args: (pangenome: Pangenome, org_name: str, contig_dict: dict, circular_contigs: dict, link: bool)

    :return: Nothing function not called yet
    """
    return read_organism(*args)


def read_organism(pangenome: Pangenome, org_name: str, contig_dict: dict, circular_contigs: dict, link: bool = False):
    """
    Read information from pangenome to assign to organism object

    :param pangenome: Input     pangenome
    :param org_name: Name of the organism
    :param contig_dict: Dictionary with all contig and associate genes
    :param circular_contigs: Dictionary of contigs
    :param link: get the gene object if the genes are clustered
    """
    org = Organism(org_name)
    gene, gene_type = (None, None)
    for contigName, geneList in contig_dict.items():
        contig = org.get_contig(contigName, is_circular=circular_contigs[contigName])
        for row in geneList:
            if link:  # if the gene families are already computed/loaded the gene exists.
                gene = pangenome.get_gene(row["ID"].decode())
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
            if isinstance(gene, Gene):
                gene.fill_annotations(start=row["start"], stop=row["stop"], strand=row["strand"].decode(),
                                      gene_type=row["type"].decode(), name=row["name"].decode(), position=row['position'],
                                      genetic_code=row["genetic_code"], product=row["product"].decode(),
                                      local_identifier=local)
            else:
                gene.fill_annotations(start=row["start"], stop=row["stop"], strand=row["strand"].decode(),
                                      gene_type=row["type"].decode(), name=row["name"].decode(),
                                      product=row["product"].decode(), local_identifier=local)
            gene.is_fragment = row["is_fragment"]
            gene.fill_parents(org, contig)
            if gene_type == "CDS":
                contig.add_gene(gene)
            elif "RNA" in gene_type:
                contig.add_rna(gene)
            else:
                raise Exception(f"A strange type '{gene_type}', which we do not know what to do with, was met.")
    pangenome.add_organism(org)


def read_graph(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read information about graph in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without graph information
    :param h5f: Pangenome HDF5 file with graph information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.edges

    if not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] or \
            not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the graph "
                        "if the annotations and the gene families have not been loaded.")
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="contig adjacency", disable=disable_bar):
        source = pangenome.get_gene(row["geneSource"].decode())
        target = pangenome.get_gene(row["geneTarget"].decode())
        pangenome.add_edge(source, target)
    pangenome.status["neighborsGraph"] = "Loaded"


def read_gene_families(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read gene families in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.geneFamilies

    link = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] else False

    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene family", disable=disable_bar):
        fam = pangenome.add_gene_family(row["geneFam"].decode())
        if link:  # linking if we have loaded the annotations
            gene_obj = pangenome.get_gene(row["gene"].decode())
        else:  # else, no
            gene_obj = Gene(row["gene"].decode())
        fam.add_gene(gene_obj)
    pangenome.status["genesClustered"] = "Loaded"


def read_gene_families_info(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read information about gene families in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families information
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.geneFamiliesInfo

    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene family", disable=disable_bar):
        fam = pangenome.add_gene_family(row["name"].decode())
        fam.add_partition(row["partition"].decode())
        fam.add_sequence(row["protein"].decode())

    if h5f.root.status._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "Loaded"
    if h5f.root.status._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "Loaded"


def read_gene_sequences(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read gene sequences in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene sequence associate to gene
    :param h5f: Pangenome HDF5 file with gene sequence associate to gene
    :param disable_bar: Disable the progress bar
    """
    if not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the pangenome gene dna sequences "
                        "if the annotations have not been loaded.")
    table = h5f.root.geneSequences

    seqid2seq = read_sequences(h5f)
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene", disable=disable_bar):
        gene = pangenome.get_gene(row['gene'].decode())
        gene.add_dna(seqid2seq[row['seqid'].decode()])
    pangenome.status["geneSequences"] = "Loaded"


def read_rgp(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read region of genomic plasticty in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without RGP
    :param h5f: Pangenome HDF5 file with RGP computed
    :param disable_bar: Disable the progress bar
    """
    if not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] or \
            not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the RGP "
                        "if the annotations and the gene families have not been loaded.")
    table = h5f.root.RGP

    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="region", disable=disable_bar):
        region = pangenome.get_region(row["RGP"].decode())
        region.append(pangenome.get_gene(row["gene"].decode()))
    # order the genes properly in the regions
    for region in pangenome.regions:
        region.genes = sorted(region.genes, key=lambda x: x.position)  # order the same way as on the contig
    pangenome.status["predictedRGP"] = "Loaded"


def read_spots(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read hotspot in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without spot
    :param h5f: Pangenome HDF5 file with spot computed
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.spots
    spots = {}
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="spot", disable=disable_bar):
        curr_spot = spots.get(row["spot"])
        if curr_spot is None:
            curr_spot = Spot(row["spot"])
            spots[row["spot"]] = curr_spot
        curr_spot.add_region(pangenome.get_region(row["RGP"].decode()))
        curr_spot.spot_2_families()
    pangenome.add_spots(spots.values())
    pangenome.status["spots"] = "Loaded"


def read_modules(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read modules in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without modules
    :param h5f: Pangenome HDF5 file with modules computed
    :param disable_bar: Disable the progress bar
    """
    if not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the modules if the gene families have not been loaded.")
    table = h5f.root.modules
    modules = {}  # id2mod
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="module", disable=disable_bar):
        curr_module = modules.get(row['module'])
        if curr_module is None:
            curr_module = Module(row['module'])
            modules[row["module"]] = curr_module
        curr_module.add_family(pangenome.get_gene_family(row['geneFam'].decode()))
    pangenome.add_modules(modules.values())
    pangenome.status["modules"] = "Loaded"


def read_annotation(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read annotation in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without annotation
    :param h5f: Pangenome HDF5 file with annotation
    :param disable_bar: Disable the progress bar
    """
    annotations = h5f.root.annotations

    table = annotations.genes
    pangenome_dict = {}
    circular_contigs = {}
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene", disable=disable_bar):
        decode_org = row["organism"].decode()
        try:
            # new gene, seen contig, seen org
            pangenome_dict[decode_org][row["contig"]["name"].decode()].append(row["gene"])
        except KeyError:
            try:
                # new contig, seen org
                pangenome_dict[decode_org][row["contig"]["name"].decode()] = [row["gene"]]
                circular_contigs[decode_org][row["contig"]["name"].decode()] = row["contig"]["is_circular"]
            except KeyError:
                # new org
                pangenome_dict[sys.intern(decode_org)] = {row["contig"]["name"].decode(): [row["gene"]]}
                circular_contigs[decode_org] = {row["contig"]["name"].decode(): row["contig"]["is_circular"]}

    link = True if pangenome.status["genesClustered"] in ["Computed", "Loaded"] else False

    for orgName, contigDict in tqdm(pangenome_dict.items(), total=len(pangenome_dict),
                                    unit="organism", disable=disable_bar):
        # TODO read organism in multiprocessing
        read_organism(pangenome, orgName, contigDict, circular_contigs[orgName], link)
    pangenome.status["genomesAnnotated"] = "Loaded"


def read_info(h5f: tables.File):
    """
    Read the pangenome content

    :param h5f: Pangenome HDF5 file
    """
    if "/info" in h5f:
        info_group = h5f.root.info

        print(f"Genes: {info_group._v_attrs['numberOfGenes']}")
        if "numberOfOrganisms" in info_group._v_attrs._f_list():
            print(f"Organisms: {info_group._v_attrs['numberOfOrganisms']}")
        if "numberOfClusters" in info_group._v_attrs._f_list():
            print(f"Families: {info_group._v_attrs['numberOfClusters']}")
        if "numberOfEdges" in info_group._v_attrs._f_list():
            print(f"Edges: {info_group._v_attrs['numberOfEdges']}")
        if 'numberOfCloud' in info_group._v_attrs._f_list():  # then all the others are there
            print(
                f"Persistent ({', '.join([key + ':' + str(round(val, 2)) for key, val in info_group._v_attrs['persistentStats'].items()])} ): "
                f"{info_group._v_attrs['numberOfPersistent']}")
            print(
                f"Shell ( {', '.join([key + ':' + str(round(val, 2)) for key, val in info_group._v_attrs['shellStats'].items()])} ): "
                f"{info_group._v_attrs['numberOfShell']}")
            print(
                f"Cloud ( {', '.join([key + ':' + str(round(val, 2)) for key, val in info_group._v_attrs['cloudStats'].items()])} ): "
                f"{info_group._v_attrs['numberOfCloud']}")
            print(f"Number of partitions: {info_group._v_attrs['numberOfPartitions']}")
            if info_group._v_attrs['numberOfPartitions'] != 3:
                for key, val in info_group._v_attrs['numberOfSubpartitions'].items():
                    print(f"Shell {key} : {val}")
        if 'genome_fluidity' in info_group._v_attrs._f_list():
            out = "Genomes fluidity: " + \
                  ", ".join(f"{subset}={round(value, 3)}" for subset, value in
                            info_group._v_attrs['genomes_fluidity'].items())
            print(out)
        if 'family_fluidity' in info_group._v_attrs._f_list():
            out = "Families fluidity: " + \
                  ", ".join(f"{subset}={round(value, 3)}" for subset, value in
                            info_group._v_attrs['families_fluidity'].items())
            print(out)
        if 'numberOfRGP' in info_group._v_attrs._f_list():
            print(f"RGPs: {info_group._v_attrs['numberOfRGP']}")
        if 'numberOfSpots' in info_group._v_attrs._f_list():
            print(f"Spots: {info_group._v_attrs['numberOfSpots']}")
        if 'numberOfModules' in info_group._v_attrs._f_list():
            if all(x in info_group._v_attrs._f_list() for x in ['CloudSpecInModules', 'ShellSpecInModules',
                                                                'numberOfFamiliesInModules']):
                read_modules_info(h5f)
            else:
                print(f"Modules: {info_group._v_attrs['numberOfModules']}")
                print(f"Families in Modules: {info_group._v_attrs['numberOfFamiliesInModules']}")


def read_modules_info(h5f: tables.File):
    """
    Read modules information in pangenome hdf5 file

    :param h5f: Pangenome HDF5 file with RGP computed
    """
    if "/info" in h5f:
        info_group = h5f.root.info
        if all(x in info_group._v_attrs._f_list() for x in ['CloudSpecInModules', 'PersistentSpecInModules',
                                                            'ShellSpecInModules', 'numberOfFamiliesInModules',
                                                            'StatOfFamiliesInModules']):
            print(f"Modules: {info_group._v_attrs['numberOfModules']}")
            print(f"Number of Families in Modules: {info_group._v_attrs['numberOfFamiliesInModules']}")
            print(f"\tPercent of Families: persistent {info_group._v_attrs['PersistentSpecInModules']['percent']},"
                  f"shell {info_group._v_attrs['ShellSpecInModules']['percent']},"
                  f"cloud {info_group._v_attrs['CloudSpecInModules']['percent']}")
            print(f"Number of Families per Modules: "
                  f"min: {info_group._v_attrs['StatOfFamiliesInModules']['min']}, "
                  f"max: {info_group._v_attrs['StatOfFamiliesInModules']['max']}, "
                  f"sd: {info_group._v_attrs['StatOfFamiliesInModules']['sd']}, "
                  f"mean: {info_group._v_attrs['StatOfFamiliesInModules']['mean']}")


def read_parameters(h5f: tables.File):
    """
    Read pangenome parameters

    :param h5f: Pangenome HDF5 file
    """
    if "/info" in h5f:
        info_group = h5f.root.info
        if "parameters" in info_group._v_attrs._f_list():
            for key, dic in info_group._v_attrs["parameters"].items():
                print(f"{key}")
                for key2, val in dic.items():
                    print(f"    {key2} : {val}")


def read_pangenome(pangenome, annotation: bool = False, gene_families: bool = False, graph: bool = False,
                   rgp: bool = False, spots: bool = False, gene_sequences: bool = False, modules: bool = False,
                   disable_bar: bool = False):
    """
    Reads a previously written pan, with all of its parts, depending on what is asked,
    with regard to what is filled in the 'status' field of the hdf5 file.

    :param pangenome: Pangenome object without some information
    :param annotation: get annotation
    :param gene_families: get gene families
    :param graph: get graph
    :param rgp: get RGP
    :param spots: get hotspot
    :param gene_sequences: get gene sequences
    :param modules: get modules
    :param disable_bar: Allow to disable the progress bar
    """
    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")
    if annotation:
        if h5f.root.status._v_attrs.genomesAnnotated:
            logging.getLogger().info("Reading pangenome annotations...")
            read_annotation(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' has not been annotated, or has been improperly filled")
    if gene_sequences:
        if h5f.root.status._v_attrs.geneSequences and h5f.root.status._v_attrs.sequences:
            logging.getLogger().info("Reading pangenome gene dna sequences...")
            read_gene_sequences(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have gene sequences, "
                            f"or has been improperly filled")

    if gene_families:
        if h5f.root.status._v_attrs.genesClustered:
            logging.getLogger().info("Reading pangenome gene families...")
            read_gene_families(pangenome, h5f, disable_bar=disable_bar)
            read_gene_families_info(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have gene families, or has been improperly filled")
    if graph:
        if h5f.root.status._v_attrs.NeighborsGraph:
            logging.getLogger().info("Reading the neighbors graph edges...")
            read_graph(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have graph information, "
                            f"or has been improperly filled")
    if rgp:
        if h5f.root.status._v_attrs.predictedRGP:
            logging.getLogger().info("Reading the RGP...")
            read_rgp(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have RGP information, "
                            f"or has been improperly filled")
    if spots:
        if h5f.root.status._v_attrs.spots:
            logging.getLogger().info("Reading the spots...")
            read_spots(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have spots information, "
                            f"or has been improperly filled")
    if modules:
        if h5f.root.status._v_attrs.modules:
            logging.getLogger().info("Reading the modules...")
            read_modules(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have modules information, "
                            f"or has been improperly filled")
    h5f.close()


def check_pangenome_info(pangenome, need_annotations: bool = False, need_families: bool = False,
                         need_graph: bool = False, need_partitions: bool = False, need_rgp: bool = False,
                         need_spots: bool = False, need_gene_sequences: bool = False, need_modules: bool = False,
                         disable_bar: bool = False):
    """
    Defines what needs to be read depending on what is needed, and automatically checks if the required elements
    have been computed with regard to the `pangenome.status`

    :param pangenome: Pangenome object without some information
    :param need_annotations: get annotation
    :param need_families: get gene families
    :param need_graph: get graph
    :param need_partitions: get partition
    :param need_rgp: get RGP
    :param need_spots: get hotspot
    :param need_gene_sequences: get gene sequences
    :param need_modules: get modules
    :param disable_bar: Allow to disable the progress bar
    """
    annotation = False
    gene_families = False
    graph = False
    rgp = False
    spots = False
    gene_sequences = False
    modules = False

    # TODO Automate call if one need another
    if need_annotations:
        if pangenome.status["genomesAnnotated"] == "inFile":
            annotation = True
        elif not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome has no genes. See the 'annotate' subcommand.")
    if need_families:
        if pangenome.status["genesClustered"] == "inFile":
            gene_families = True
        elif not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome has no gene families. See the 'cluster' subcommand.")
    if need_graph:
        if pangenome.status["neighborsGraph"] == "inFile":
            graph = True
        elif not pangenome.status["neighborsGraph"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome does not have a graph (no edges). See the 'graph' subcommand.")
    if need_partitions and pangenome.status["partitioned"] not in ["Computed", "Loaded", "inFile"]:
        raise Exception("Your pangenome has not been partitioned. See the 'partition' subcommand")
    if need_rgp:
        if pangenome.status["predictedRGP"] == "inFile":
            rgp = True
        elif not pangenome.status["predictedRGP"] in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome  regions of genomic plasticity have not been predicted. See the 'rgp' subcommand")
    if need_spots:
        if pangenome.status["spots"] == "inFile":
            spots = True
        elif not pangenome.status["spots"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome spots of insertion have not been predicted. See the 'spot' subcommand")
    if need_gene_sequences:
        if pangenome.status["geneSequences"] == "inFile":
            gene_sequences = True
        elif not pangenome.status["geneSequences"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome does not include gene sequences. "
                            "This is possible only if you provided your own cluster file with the 'cluster' subcommand")

    if need_modules:
        if pangenome.status["modules"] == "inFile":
            modules = True
        elif not pangenome.status["modules"] in ["Computed", "Loaded"]:
            raise Exception("Your pangenome modules have not been predicted. See the 'module' subcommand")

    if annotation or gene_families or graph or rgp or spots or gene_sequences or modules:
        # if anything is true, else we need nothing.
        read_pangenome(pangenome, annotation=annotation, gene_families=gene_families, graph=graph, rgp=rgp, spots=spots,
                       gene_sequences=gene_sequences, modules=modules, disable_bar=disable_bar)
