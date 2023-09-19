#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import sys
from pathlib import Path
from typing import TextIO, List, Dict, Tuple

# installed libraries
from tables import Table
from tqdm import tqdm
import tables

# local libraries
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot, Module
from ppanggolin.metadata import Metadata



class Genedata:
    """
    This is a general class storing unique gene-related data to be written in a specific
    genedata table
    """

    def __init__(self, start: int, stop: int, strand: str, gene_type: str, position: int, name: str, product: str,
                 genetic_code: int):
        """Constructor method

        :param start: Gene start position
        :param stop: Gene stop position
        :param strand: Associated strand
        :param gene_type: Gene type
        :param position: Position of the gene on its contig
        :param name: Name of the feature
        :param product: Associated product
        :param genetic_code: associated genetic code, if any
        """
        self.start = start
        self.stop = stop
        self.strand = strand
        self.gene_type = gene_type
        self.position = position
        self.name = name
        self.product = product
        self.genetic_code = genetic_code

    def __eq__(self, other):
        return self.start == other.start \
            and self.stop == other.stop \
            and self.strand == other.strand \
            and self.gene_type == other.gene_type \
            and self.position == other.position \
            and self.name == other.name \
            and self.product == other.product \
            and self.genetic_code == other.genetic_code

    def __hash__(self):
        return hash((self.start, self.stop, self.strand, self.gene_type, self.position,
                     self.name, self.product, self.genetic_code))


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


# TODO Remove this function
def fix_partitioned(pangenome_file: str):
    """
    Fixes pangenomes with the 'partitionned' typo.

    :param pangenome_file: path to the pangenome file
    """
    h5f = tables.open_file(pangenome_file, "a")
    status_group = h5f.root.status
    if 'Partitionned' in status_group._v_attrs._f_list():
        # if Partitionned is still in use, fix it
        status_group = h5f.root.status
        if status_group._v_attrs.Partitionned:
            status_group._v_attrs.Partitioned = True
        else:
            status_group._v_attrs.Partitioned = False
        del status_group._v_attrs.Partitionned
    h5f.close()

def get_status(pangenome: Pangenome, pangenome_file: Path):
    """
    Checks which elements are already present in the file.

    :param pangenome: Blank pangenome
    :param pangenome_file: path to the pangenome file
    """
    fix_partitioned(pangenome_file)
    h5f = tables.open_file(pangenome_file, "r")
    logging.getLogger("PPanGGOLiN").info("Getting the current pangenome status")
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

    if status_group._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "inFile"

    if hasattr(status_group._v_attrs, "predictedRGP") and status_group._v_attrs.predictedRGP:
        pangenome.status["predictedRGP"] = "inFile"

    if hasattr(status_group._v_attrs, "spots") and status_group._v_attrs.spots:
        pangenome.status["spots"] = "inFile"

    if hasattr(status_group._v_attrs, "modules") and status_group._v_attrs.modules:
        pangenome.status["modules"] = "inFile"
        # pangenome.status["annotations_sources"] = status_group._v_attrs.annotations_sources

    if hasattr(status_group._v_attrs, "metadata") and status_group._v_attrs.metadata:
        metastatus = status_group.metastatus
        metasources = status_group.metasources
        for attr in metastatus._v_attrs._f_list():
            pangenome.status["metadata"][attr] = "inFile"
            pangenome.status["metasources"][attr] = metasources._v_attrs[attr]

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


def read_genedata(h5f: tables.File) -> Dict[int, Genedata]:
    """
    Reads the genedata table and returns a genedata_id2genedata dictionnary

    :param h5f: the hdf5 file handler

    :return: dictionnary linking genedata to the genedata identifier
    """
    table = h5f.root.annotations.genedata
    genedata_id2genedata = {}
    for row in read_chunks(table, chunk=20000):
        genedata = Genedata(start=int(row["start"]),
                            stop=int(row["stop"]),
                            strand=row["strand"].decode(),
                            gene_type=row["gene_type"].decode(),
                            position=int(row["position"]),
                            name=row["name"].decode(),
                            product=row["product"].decode(),
                            genetic_code=int(row["genetic_code"]))
        genedata_id = row["genedata_id"]
        genedata_id2genedata[genedata_id] = genedata
    return genedata_id2genedata


def read_sequences(h5f: tables.File) -> dict:
    """
    Reads the sequences table and returns a sequence id to sequence dictionnary
    :param h5f: the hdf5 file handler
    :return: dictionnary linking sequences to the seq identifier
    """
    table = h5f.root.annotations.sequences
    seqid2seq = {}
    for row in read_chunks(table, chunk=20000):
        seqid2seq[row["seqid"]] = row['dna'].decode()
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
    logging.getLogger("PPanGGOLiN").info(f"Extracting and writing CDS sequences from a {filename} file to a fasta file...")
    h5f = tables.open_file(filename, "r", driver_core_backing_store=0)
    table = h5f.root.annotations.geneSequences
    list_cds = set(list_cds) if list_cds is not None else None
    seqid2seq = read_sequences(h5f)
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene", disable=disable_bar):
        # Read the table chunk per chunk otherwise RAM dies on big pangenomes
        name_cds = row["gene"].decode()
        if row["type"] == b"CDS" and (list_cds is None or name_cds in list_cds):
            file_obj.write('>' + add + name_cds + "\n")
            file_obj.write(seqid2seq[row["seqid"]] + "\n")
    file_obj.flush()
    h5f.close()


def read_graph(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read information about graph in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without graph information
    :param h5f: Pangenome HDF5 file with graph information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.edges

    if pangenome.status["genomesAnnotated"] not in ["Computed", "Loaded"] or \
            pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
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
        try:
            fam = pangenome.get_gene_family(name=row["geneFam"].decode())
        except KeyError:
            fam = GeneFamily(family_id=pangenome.max_fam_id, name=row["geneFam"].decode())
            pangenome.add_gene_family(fam)
        if link:  # linking if we have loaded the annotations
            gene_obj = pangenome.get_gene(row["gene"].decode())
        else:  # else, no
            gene_obj = Gene(row["gene"].decode())
        fam.add(gene_obj)
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
        fam = pangenome.get_gene_family(row["name"].decode())
        fam.partition = row["partition"].decode()
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
    if pangenome.status["genomesAnnotated"] not in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the pangenome gene dna sequences "
                        "if the annotations have not been loaded.")
    table = h5f.root.annotations.geneSequences

    seqid2seq = read_sequences(h5f)
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene", disable=disable_bar):
        gene = pangenome.get_gene(row['gene'].decode())
        gene.add_sequence(seqid2seq[row['seqid']])
    pangenome.status["geneSequences"] = "Loaded"


def read_rgp(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read region of genomic plasticty in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without RGP
    :param h5f: Pangenome HDF5 file with RGP computed
    :param disable_bar: Disable the progress bar
    """
    if pangenome.status["genomesAnnotated"] not in ["Computed", "Loaded"] or \
            pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the RGP "
                        "if the annotations and the gene families have not been loaded.")
    table = h5f.root.RGP

    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="region", disable=disable_bar):
        try:
            region = pangenome.get_region(row["RGP"].decode())
        except KeyError:
            region = Region(row["RGP"].decode())
            pangenome.add_region(region)
        gene = pangenome.get_gene(row["gene"].decode())
        region.add(gene)
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
        curr_spot = spots.get(int(row["spot"]))
        if curr_spot is None:
            curr_spot = Spot(int(row["spot"]))
            spots[row["spot"]] = curr_spot
        region = pangenome.get_region(row["RGP"].decode())
        curr_spot.add(region)
        curr_spot.spot_2_families()
    for spot in spots.values():
        pangenome.add_spot(spot)
    pangenome.status["spots"] = "Loaded"


def read_modules(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read modules in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without modules
    :param h5f: Pangenome HDF5 file with modules computed
    :param disable_bar: Disable the progress bar
    """
    if pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the modules if the gene families have not been loaded.")
    table = h5f.root.modules
    modules = {}  # id2mod
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="module", disable=disable_bar):
        curr_module = modules.get(int(row['module']))
        if curr_module is None:
            curr_module = Module(int(row['module']))
            modules[row["module"]] = curr_module
        family = pangenome.get_gene_family(row['geneFam'].decode())
        curr_module.add(family)
    for module in modules.values():
        pangenome.add_module(module)
    pangenome.status["modules"] = "Loaded"

def read_organisms(pangenome: Pangenome, table: tables.Table, chunk_size: int = 20000,
                   disable_bar: bool = False):
    """Read organism table in pangenome file to add them to the pangenome object

    :param pangenome: Pangenome object
    :param table: Organism table
    :param chunk_size: Size of the chunck reading
    :param disable_bar: Disable progress bar
    """
    contig2organism = {}
    for row in tqdm(read_chunks(table, chunk=chunk_size), total=table.nrows, unit="genome", disable=disable_bar):
        organism = Organism(row["name"].decode())
        pangenome.add_organism(organism)


def read_contigs(pangenome: Pangenome, table: tables.Table, chunk_size: int = 20000,
                   disable_bar: bool = False):
    """Read contig table in pangenome file to add them to the pangenome object

    :param pangenome: Pangenome object
    :param table: Contig table
    :param chunk_size: Size of the chunck reading
    :param disable_bar: Disable progress bar
    """
    for row in tqdm(read_chunks(table, chunk=chunk_size), total=table.nrows, unit="contig", disable=disable_bar):
        contig = Contig(name=row["name"].decode())
        contig.is_circular = row["is_circular"]
        contig.length = int(row["length"])
        try:
            organism = pangenome.get_organism(row["organism"].decode())
        except KeyError:
            pass
        else:
            organism.add(contig)

def read_genes(pangenome: Pangenome, table: tables.Table, genedata_dict: Dict[int, Genedata],
               link: bool = True, chunk_size: int = 20000, disable_bar: bool = False):
    """Read genes in pangenome file to add them to the pangenome object

    :param pangenome: Pangenome object
    :param table: Genes table
    :param genedata_dict: Dictionary to link genedata with gene
    :param link: Allow to link gene to organism and contig
    :param chunk_size: Size of the chunck reading
    :param disable_bar: Disable progress bar
    """
    for row in tqdm(read_chunks(table, chunk=chunk_size), total=table.nrows, unit="gene", disable=disable_bar):
        gene = Gene(row["ID"].decode())
        genedata = genedata_dict[row["genedata_id"]]
        try:
            local = row["local"].decode()
        except ValueError:
            local = ""
        gene.fill_annotations(start=genedata.start, stop=genedata.stop, strand=genedata.strand,
                              gene_type=genedata.gene_type, name=genedata.name, position=genedata.position,
                              genetic_code=genedata.genetic_code, product=genedata.product, local_identifier=local)
        gene.is_fragment = row["is_fragment"]
        if link:
            contig = pangenome.get_contig(row["contig"].decode())
            gene.fill_parents(contig.organism, contig)
            contig.add(gene)


def read_rnas(pangenome: Pangenome, table: tables.Table, genedata_dict: Dict[int, Genedata],
              link: bool = True, chunk_size: int = 20000, disable_bar: bool = False):
    """Read RNAs in pangenome file to add them to the pangenome object

    :param pangenome: Pangenome object
    :param table: RNAs table
    :param genedata_dict: Dictionary to link genedata with gene
    :param link: Allow to link gene to organism and contig
    :param chunk_size: Size of the chunck reading
    :param disable_bar: Disable progress bar
    """
    for row in tqdm(read_chunks(table, chunk=chunk_size), total=table.nrows, unit="gene", disable=disable_bar):
        rna = RNA(row["ID"].decode())
        genedata = genedata_dict[row["genedata_id"]]
        rna.fill_annotations(start=genedata.start, stop=genedata.stop, strand=genedata.strand,
                             gene_type=genedata.gene_type, name=genedata.name,
                             product=genedata.product)
        if link:
            contig = pangenome.get_contig(row["contig"].decode())
            rna.fill_parents(contig.organism, contig)
            contig.add_rna(rna)


def read_annotation(pangenome: Pangenome, h5f: tables.File, load_organisms: bool = True, load_contigs: bool = True,
                    load_genes: bool = True, load_rnas: bool = True, chunk_size: int = 20000, disable_bar: bool = False):
    """
    Read annotation in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without annotation
    :param h5f: Pangenome HDF5 file with annotation
    :param disable_bar: Disable the progress bar
    """
    annotations = h5f.root.annotations
    genedata_dict = None
    if load_organisms:
        read_organisms(pangenome, annotations.genomes, disable_bar=disable_bar)

    if load_contigs:
        read_contigs(pangenome, annotations.contigs, disable_bar=disable_bar)

    if load_genes:
        genedata_dict = read_genedata(h5f)
        read_genes(pangenome, annotations.genes, genedata_dict,
                   all([load_organisms, load_contigs]), disable_bar=disable_bar)
    if load_rnas:
        read_rnas(pangenome, annotations.RNAs, read_genedata(h5f) if genedata_dict is None else genedata_dict,
                   all([load_organisms, load_contigs]), disable_bar=disable_bar)
    pangenome.status["genomesAnnotated"] = "Loaded"



def read_info(h5f: tables.File):
    """
    Read the pangenome content

    :param h5f: Pangenome HDF5 file
    """
    if "/info" in h5f:
        info_group = h5f.root.info
        print("Content: ")
        print(f"\t- Genes: {info_group._v_attrs['numberOfGenes']}")
        if "numberOfOrganisms" in info_group._v_attrs._f_list():
            print(f"\t- Organisms: {info_group._v_attrs['numberOfOrganisms']}")
        if "numberOfClusters" in info_group._v_attrs._f_list():
            print(f"\t- Families: {info_group._v_attrs['numberOfClusters']}")
        if "numberOfEdges" in info_group._v_attrs._f_list():
            print(f"\t- Edges: {info_group._v_attrs['numberOfEdges']}")
        if 'numberOfCloud' in info_group._v_attrs._f_list():  # then all the others are there
            print(f"\t- Persistent: \n"
                  f"\t\t- count : {info_group._v_attrs['numberOfPersistent']}")
            for key, val in info_group._v_attrs['persistentStats'].items():
                print(f"\t\t- {key}: {str(round(val, 2))}")
            print(f"\t- Shell: \n"
                  f"\t\t- count : {info_group._v_attrs['numberOfShell']}")
            for key, val in info_group._v_attrs['shellStats'].items():
                print(f"\t\t- {key}: {str(round(val, 2))}")
            print(f"\t- Cloud: \n"
                  f"\t\t- count : {info_group._v_attrs['numberOfCloud']}")
            for key, val in info_group._v_attrs['cloudStats'].items():
                print(f"\t\t- {key}: {str(round(val, 2))}")
            print(f"\t- Number of partitions: {info_group._v_attrs['numberOfPartitions']}")
            if info_group._v_attrs['numberOfPartitions'] != 3:
                for key, val in info_group._v_attrs['numberOfSubpartitions'].items():
                    print(f"\t\t- Shell {key} : {val}")
        if 'genomes_fluidity' in info_group._v_attrs._f_list():
            print("\t- Genomes fluidity: ")
            for subset, value in info_group._v_attrs['genomes_fluidity'].items():
                print(f"\t\t- {subset}: {round(value, 3)}")
        if 'family_fluidity' in info_group._v_attrs._f_list():
            out = "\t- Families fluidity: " + \
                  ", ".join(f"{subset}={round(value, 3)}" for subset, value in
                            info_group._v_attrs['families_fluidity'].items())
            print(out)
        if 'numberOfRGP' in info_group._v_attrs._f_list():
            print(f"\t- RGPs: {info_group._v_attrs['numberOfRGP']}")
        if 'numberOfSpots' in info_group._v_attrs._f_list():
            print(f"\t- Spots: {info_group._v_attrs['numberOfSpots']}")
        if 'numberOfModules' in info_group._v_attrs._f_list():
            if all(x in info_group._v_attrs._f_list() for x in ['CloudSpecInModules', 'ShellSpecInModules',
                                                                'numberOfFamiliesInModules']):
                read_modules_info(h5f)
            else:
                print(f"\t- Modules: {info_group._v_attrs['numberOfModules']}")
                print(f"\t- Families in Modules: {info_group._v_attrs['numberOfFamiliesInModules']}")


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
            print(f"\t- Modules: {info_group._v_attrs['numberOfModules']}")
            print(f"\t\t- Families in Modules: {info_group._v_attrs['numberOfFamiliesInModules']}")
            print(f"\t\t- Percent of Families: \n"
                  f"\t\t\t- persistent: {info_group._v_attrs['PersistentSpecInModules']['percent']}\n"
                  f"\t\t\t- shell {info_group._v_attrs['ShellSpecInModules']['percent']}\n"
                  f"\t\t\t- cloud {info_group._v_attrs['CloudSpecInModules']['percent']}")
            print(f"\t\t- Number of Families per Modules:\n"
                  f"\t\t\t- min: {info_group._v_attrs['StatOfFamiliesInModules']['min']}\n"
                  f"\t\t\t- max: {info_group._v_attrs['StatOfFamiliesInModules']['max']}\n"
                  f"\t\t\t- sd: {info_group._v_attrs['StatOfFamiliesInModules']['sd']}\n"
                  f"\t\t\t- mean: {info_group._v_attrs['StatOfFamiliesInModules']['mean']}")


def read_metadata(pangenome: Pangenome, h5f: tables.File, metatype: str,
                  sources: List[str] = None, disable_bar: bool = False):
    """Read metadata to add them to the pangenome object

    :param pangenome: Pangenome object
    :param h5f: Pangenome file
    :param metatype: Object type to associate metadata
    :param sources: Source name of metadata
    :param disable_bar: Disable progress bar
    """
    metadata_group = h5f.root.metadata._f_get_child(metatype)
    for source in sources:
        source_table = metadata_group._f_get_child(source)
        for row in tqdm(read_chunks(source_table), total=source_table.nrows, unit='metadata', disable=disable_bar):
            meta_dict = {'source': source}
            if "ID" in row.dtype.names:
                identifier = row["ID"].decode() if isinstance(row["ID"], bytes) else row["ID"]
            else:
                identifier = row["name"].decode()
            if metatype == "families":
                element = pangenome.get_gene_family(identifier)
            elif metatype == "genomes":
                element = pangenome.get_organism(identifier)
            elif metatype == "genes":
                element = pangenome.get_gene(identifier)
            elif metatype == "RGPs":
                element = pangenome.get_region(identifier)
            elif metatype == "spots":
                element = pangenome.get_spot(identifier)
            else:  # metatype == "modules":
                element = pangenome.get_module(identifier)
            for field in row.dtype.names:
                if field not in ["ID", "name"]:
                    meta_dict[field] = row[field].decode() if isinstance(row[field], bytes) else row[field]
            meta = Metadata(**meta_dict)
            element.add_metadata(source=source, metadata=meta)
    pangenome.status["metadata"][metatype] = "Loaded"

def read_parameters(h5f: tables.File):
    """
    Read pangenome parameters

    :param h5f: Pangenome HDF5 file
    """
    if "/info" in h5f:
        info_group = h5f.root.info
        if "parameters" in info_group._v_attrs._f_list():
            print("Parameters: ")
            for key, dic in info_group._v_attrs["parameters"].items():
                print(f"\t- {key}")
                for key2, val in dic.items():
                    print(f"\t\t- {key2} : {val}")


def read_pangenome(pangenome, annotation: bool = False, gene_families: bool = False, graph: bool = False,
                   rgp: bool = False, spots: bool = False, gene_sequences: bool = False, modules: bool = False,
                   metadata: bool = False, metatype: str = None, sources: List[str] = None,
                   disable_bar: bool = False):
    """
    Reads a previously written pangenome, with all of its parts, depending on what is asked,
    with regard to what is filled in the 'status' field of the hdf5 file.

    :param pangenome: Pangenome object without some information
    :param annotation: get annotation
    :param gene_families: get gene families
    :param graph: get graph
    :param rgp: get RGP
    :param spots: get hotspot
    :param gene_sequences: get gene sequences
    :param modules: get modules
    :param metadata: get metadata
    :param metatype: metatype of the metadata to get
    :param sources: sources of the metadata to get (None means all sources)
    :param disable_bar: Allow to disable the progress bar
    """
    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")

    fix_partitioned(pangenome.file)

    h5f = tables.open_file(filename, "r")

    if annotation:  # I place annotation here, to link gene to gene families if organism are not loaded
        if h5f.root.status._v_attrs.genomesAnnotated:
            logging.getLogger("PPanGGOLiN").info("Reading pangenome annotations...")
            read_annotation(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' has not been annotated, or has been improperly filled")

    if gene_sequences:
        if h5f.root.status._v_attrs.geneSequences:
            logging.getLogger("PPanGGOLiN").info("Reading pangenome gene dna sequences...")
            read_gene_sequences(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have gene sequences, "
                            f"or has been improperly filled")

    if gene_families:
        if h5f.root.status._v_attrs.genesClustered:
            logging.getLogger("PPanGGOLiN").info("Reading pangenome gene families...")
            read_gene_families(pangenome, h5f, disable_bar=disable_bar)
            read_gene_families_info(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have gene families, or has been improperly filled")

    if graph:
        if h5f.root.status._v_attrs.NeighborsGraph:
            logging.getLogger("PPanGGOLiN").info("Reading the neighbors graph edges...")
            read_graph(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have graph information, "
                            f"or has been improperly filled")

    if rgp:
        if h5f.root.status._v_attrs.predictedRGP:
            logging.getLogger("PPanGGOLiN").info("Reading the RGP...")
            read_rgp(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have RGP information, "
                            f"or has been improperly filled")

    if spots:
        if h5f.root.status._v_attrs.spots:
            logging.getLogger("PPanGGOLiN").info("Reading the spots...")
            read_spots(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have spots information, "
                            f"or has been improperly filled")

    if modules:
        if h5f.root.status._v_attrs.modules:
            logging.getLogger("PPanGGOLiN").info("Reading the modules...")
            read_modules(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have modules information, "
                            f"or has been improperly filled")

    if metadata:
        assert metatype is not None
        if sources is None:
            sources = pangenome.status["metasources"][metatype]
        if h5f.root.status._v_attrs.metadata:
            metastatus = h5f.root.status._f_get_child("metastatus")
            metasources = h5f.root.status._f_get_child("metasources")
            if metastatus._v_attrs[metatype] and all([True if source in metasources._v_attrs[metatype] else False for source in sources]):
                logging.getLogger().info(f"Reading the {metatype} metadata from sources {sources}...")
                read_metadata(pangenome, h5f, metatype, sources, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' does not have modules information, "
                            f"or has been improperly filled")
    h5f.close()

def check_pangenome_info(pangenome, need_annotations: bool = False, need_families: bool = False,
                         need_graph: bool = False, need_partitions: bool = False, need_rgp: bool = False,
                         need_spots: bool = False, need_gene_sequences: bool = False, need_modules: bool = False,
                         need_metadata: bool = False, metatype: str = None, sources: List[str] = None,
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
    :param need_metadata: get metadata
    :param metatype: metatype of the metadata to get
    :param sources: sources of the metadata to get (None means all sources)
    :param disable_bar: Allow to disable the progress bar
    """
    annotation = False
    gene_families = False
    graph = False
    rgp = False
    spots = False
    gene_sequences = False
    modules = False
    metadata = False

    # TODO Automate call if one need another
    if need_annotations:
        if pangenome.status["genomesAnnotated"] == "inFile":
            annotation = True
        elif pangenome.status["genomesAnnotated"] not in ["Computed", "Loaded"]:
            raise Exception("Your pangenome has no genes. See the 'annotate' subcommand.")
    if need_families:
        if pangenome.status["genesClustered"] == "inFile":
            gene_families = True
        elif pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
            raise Exception("Your pangenome has no gene families. See the 'cluster' subcommand.")
    if need_graph:
        if pangenome.status["neighborsGraph"] == "inFile":
            graph = True
        elif pangenome.status["neighborsGraph"] not in ["Computed", "Loaded"]:
            raise Exception("Your pangenome does not have a graph (no edges). See the 'graph' subcommand.")
    if need_partitions and pangenome.status["partitioned"] not in ["Computed", "Loaded", "inFile"]:
        raise Exception("Your pangenome has not been partitioned. See the 'partition' subcommand")
    if need_rgp:
        if pangenome.status["predictedRGP"] == "inFile":
            rgp = True
        elif pangenome.status["predictedRGP"] not in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome  regions of genomic plasticity have not been predicted. See the 'rgp' subcommand")
    if need_spots:
        if pangenome.status["spots"] == "inFile":
            spots = True
        elif pangenome.status["spots"] not in ["Computed", "Loaded"]:
            raise Exception("Your pangenome spots of insertion have not been predicted. See the 'spot' subcommand")
    if need_gene_sequences:
        if pangenome.status["geneSequences"] == "inFile":
            gene_sequences = True
        elif pangenome.status["geneSequences"] not in ["Computed", "Loaded"]:
            raise Exception("Your pangenome does not include gene sequences. "
                            "This is possible only if you provided your own cluster file with the 'cluster' subcommand")

    if need_modules:
        if pangenome.status["modules"] == "inFile":
            modules = True
        elif pangenome.status["modules"] not in ["Computed", "Loaded"]:
            raise Exception("Your pangenome modules have not been predicted. See the 'module' subcommand")

    if need_metadata:
        if pangenome.status["metadata"][metatype] == "inFile":
            if sources is not None:
                for source in sources:
                    if source in pangenome.status["metasources"][metatype]:
                        metadata = True
                    else:
                        raise Exception(f"There is no metadata assign to {metatype} for source : {source} in your pangenome.")
            else:
                metadata = True
        elif not pangenome.status["metastatus"][metatype] in ["Computed", "Loaded"]:
            raise Exception(f"Your pangenome don't have any metadata for {metatype}. See the 'metadata' subcommand")

    if any([annotation, gene_families, graph, rgp, spots, gene_sequences, modules, metadata]):
        # if anything is true, else we need nothing.
        read_pangenome(pangenome, annotation=annotation, gene_families=gene_families,
                       graph=graph, gene_sequences=gene_sequences,
                       rgp=rgp, spots=spots, modules=modules,
                       metadata=metadata, metatype=metatype, sources=sources,
                       disable_bar=disable_bar)
