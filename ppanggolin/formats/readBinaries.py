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


def get_number_of_organisms(pangenome):
    """ standalone function to get the number of organisms in a pan"""
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


def get_status(pangenome, pangenome_file):
    """
        Checks which elements are already present in the file.
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
        if status_group._v_attrs.Partitionned:
            status_group._v_attrs.Partitioned = True
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


def read_chunks(table, column=None, chunk=10000):
    """
        Reading entirely the provided table (or column if specified) chunk per chunk to limit RAM usage.
    """
    for i in range(0, table.nrows, chunk):
        for row in table.read(start=i, stop=i + chunk, field=column):
            yield row


def get_gene_sequences_from_file(filename, file_obj, list_cds=None, add='', disable_bar=False):
    """
        Writes the CDS sequences of the Pangenome object to a File object that can be filtered or not by a list of CDS,
        and adds the eventual str 'add' in front of the identifiers
        Loads the sequences from a .h5 pangenome file
    """
    logging.getLogger().info("Extracting and writing CDS sequences from a .h5 pangenome file to a fasta file...")
    h5f = tables.open_file(filename, "r", driver_core_backing_store=0)
    table = h5f.root.geneSequences
    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    list_cds = set(list_cds) if list_cds is not None else None
    for row in read_chunks(table,
                           chunk=20000):  # reading the table chunk per chunk otherwise RAM dies on big pangenomes
        name_cds = row["gene"].decode()
        if row["type"] == b"CDS" and (list_cds is None or name_cds in list_cds):
            file_obj.write('>' + add + name_cds + "\n")
            file_obj.write(row["dna"].decode() + "\n")
        bar.update()
    file_obj.flush()
    bar.close()
    h5f.close()


def launch_read_organism(args):
    return read_organism(*args)


def read_organism(pangenome, org_name, contig_dict, circular_contigs, link=False):
    org = Organism(org_name)
    gene, gene_type = (None, None)
    for contigName, geneList in contig_dict.items():
        contig = org.get_or_add_contig(contigName, is_circular=circular_contigs[contigName])
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
            gene.fill_annotations(start=row["start"], stop=row["stop"], strand=row["strand"].decode(),
                                  gene_type=row["type"].decode(), name=row["name"].decode(),
                                  product=row["product"].decode(), local_identifier=local, position=row["position"],
                                  genetic_code=row["genetic_code"])
            gene.is_fragment = row["is_fragment"]
            gene.fill_parents(org, contig)
            if gene_type == "CDS":
                contig.add_gene(gene)
            elif "RNA" in gene_type:
                contig.add_rna(gene)
            else:
                raise Exception(f"A strange type '{gene_type}', which we do not know what to do with, was met.")
    pangenome.add_organism(org)


def read_graph(pangenome, h5f, disable_bar=False):
    table = h5f.root.edges

    if not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] or \
            not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the graph "
                        "if the annotations and the gene families have not been loaded.")

    bar = tqdm(range(table.nrows), unit="contig adjacency", disable=disable_bar)
    for row in read_chunks(table):
        source = pangenome.get_gene(row["geneSource"].decode())
        target = pangenome.get_gene(row["geneTarget"].decode())
        pangenome.add_edge(source, target)
        bar.update()
    bar.close()
    pangenome.status["neighborsGraph"] = "Loaded"


def read_gene_families(pangenome, h5f, disable_bar=False):
    table = h5f.root.gene_families

    link = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] else False

    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    for row in read_chunks(table):
        fam = pangenome.add_gene_family(row["geneFam"].decode())
        if link:  # linking if we have loaded the annotations
            gene_obj = pangenome.get_gene(row["gene"].decode())
        else:  # else, no
            gene_obj = Gene(row["gene"].decode())
        fam.add_gene(gene_obj)
        bar.update()
    bar.close()
    pangenome.status["genesClustered"] = "Loaded"


def read_gene_families_info(pangenome, h5f, disable_bar=False):
    table = h5f.root.geneFamiliesInfo

    bar = tqdm(range(table.nrows), unit="gene family", disable=disable_bar)
    for row in read_chunks(table):
        fam = pangenome.add_gene_family(row["name"].decode())
        fam.add_partition(row["partition"].decode())
        fam.add_sequence(row["protein"].decode())
        bar.update()
    bar.close()

    if h5f.root.status._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "Loaded"
    if h5f.root.status._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "Loaded"


def read_gene_sequences(pangenome, h5f, disable_bar=False):
    if not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the pangenome gene dna sequences "
                        "if the annotations have not been loaded.")
    table = h5f.root.geneSequences

    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    for row in read_chunks(table):
        gene = pangenome.get_gene(row['gene'].decode())
        gene.add_dna(row['dna'].decode())
        bar.update()
    bar.close()
    pangenome.status["geneSequences"] = "Loaded"


def read_rgp(pangenome, h5f, disable_bar=False):
    if not pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] or \
            not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the RGP "
                        "if the annotations and the gene families have not been loaded.")
    table = h5f.root.RGP

    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    for row in read_chunks(table):
        region = pangenome.get_or_add_region(row["RGP"].decode())
        region.append(pangenome.get_gene(row["gene"].decode()))
        bar.update()
    bar.close()
    # order the genes properly in the regions
    for region in pangenome.regions:
        region.genes = sorted(region.genes, key=lambda x: x.position)  # order the same way as on the contig
    pangenome.status["predictedRGP"] = "Loaded"


def read_spots(pangenome, h5f, disable_bar=False):
    table = h5f.root.spots
    bar = tqdm(range(table.nrows), unit="region", disable=disable_bar)
    spots = {}
    for row in read_chunks(table):
        curr_spot = spots.get(row["spot"])
        if curr_spot is None:
            curr_spot = Spot(row["spot"])
            spots[row["spot"]] = curr_spot
        curr_spot.add_region(pangenome.get_or_add_region(row["RGP"].decode()))
        curr_spot.spot_2_families()
        bar.update()
    bar.close()
    pangenome.add_spots(spots.values())
    pangenome.status["spots"] = "Loaded"


def read_modules(pangenome, h5f, disable_bar=False):
    if not pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the modules if the gene families have not been loaded.")
    table = h5f.root.modules
    bar = tqdm(range(table.nrows), unit="module", disable=disable_bar)
    modules = {}  # id2mod
    for row in read_chunks(table):
        curr_module = modules.get(row['module'])
        if curr_module is None:
            curr_module = Module(row['module'])
            modules[row["module"]] = curr_module
        curr_module.add_family(pangenome.get_gene_family(row['geneFam'].decode()))
        bar.update()
    bar.close()
    pangenome.add_modules(modules.values())
    pangenome.status["modules"] = "Loaded"


def read_annotation(pangenome, h5f, disable_bar=False):
    annotations = h5f.root.annotations

    table = annotations.genes
    bar = tqdm(range(table.nrows), unit="gene", disable=disable_bar)
    pangenome_dict = {}
    circular_contigs = {}
    for row in read_chunks(table):
        try:
            pangenome_dict[row["organism"].decode()][row["contig"]["name"].decode()].append(
                row["gene"])  # new gene, seen contig, seen org
        except KeyError:
            try:
                pangenome_dict[row["organism"].decode()][row["contig"]["name"].decode()] = [
                    row["gene"]]  # new contig, seen org
                circular_contigs[row["organism"].decode()][row["contig"]["name"].decode()] = \
                    row["contig"]["is_circular"]
            except KeyError:
                pangenome_dict[sys.intern(row["organism"].decode())] = {
                    row["contig"]["name"].decode(): [row["gene"]]}  # new org
                circular_contigs[row["organism"].decode()] = {
                    row["contig"]["name"].decode(): row["contig"]["is_circular"]}
        bar.update()
    bar.close()

    link = True if pangenome.status["genesClustered"] in ["Computed", "Loaded"] else False

    bar = tqdm(range(len(pangenome_dict)), unit="organism", disable=disable_bar)
    for orgName, contigDict in pangenome_dict.items():
        read_organism(pangenome, orgName, contigDict, circular_contigs[orgName], link)
        bar.update()
    bar.close()
    pangenome.status["genomesAnnotated"] = "Loaded"


def read_info(h5f):
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
                # readModulesInfo(h5f)


def read_modules_info(h5f):
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


def read_parameters(h5f):
    if "/info" in h5f:
        info_group = h5f.root.info
        if "parameters" in info_group._v_attrs._f_list():
            for key, dic in info_group._v_attrs["parameters"].items():
                print(f"{key}")
                for key2, val in dic.items():
                    print(f"    {key2} : {val}")


def read_pangenome(pangenome, annotation=False, gene_families=False, graph=False, rgp=False, spots=False,
                   gene_sequences=False, modules=False, disable_bar=False):
    """
        Reads a previously written pan, with all of its parts, depending on what is asked,
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
            read_annotation(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(f"The pangenome in file '{filename}' has not been annotated, or has been improperly filled")
    if gene_sequences:
        if h5f.root.status._v_attrs.geneSequences:
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
            raise Exception(
                f"The pangenome in file '{filename}' does not have graph information, or has been improperly filled")
    if rgp:
        if h5f.root.status._v_attrs.predictedRGP:
            logging.getLogger().info("Reading the RGP...")
            read_rgp(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have RGP information, or has been improperly filled")
    if spots:
        if h5f.root.status._v_attrs.spots:
            logging.getLogger().info("Reading the spots...")
            read_spots(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have spots information, or has been improperly filled")
    if modules:
        if h5f.root.status._v_attrs.modules:
            logging.getLogger().info("Reading the modules...")
            read_modules(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have modules information, or has been improperly filled")
    h5f.close()


def check_pangenome_info(pangenome, need_annotations=False, need_families=False, need_graph=False,
                         need_partitions=False, need_rgp=False, need_spots=False, need_gene_sequences=False,
                         need_modules=False, disable_bar=False):
    """
    defines what needs to be read depending on what is needed, and automatically checks if the required elements
    have been computed with regard to the `pangenome.status`
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
