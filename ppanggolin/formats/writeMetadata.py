#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from typing import List, Tuple

# installed libraries
from tqdm import tqdm
import tables

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Organism, Gene
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot, Module


def write_metadata_status(pangenome: Pangenome, h5f: tables.File, status_group: tables.Group):
    metastatus = pangenome.status["metadata"]
    metasources = pangenome.status["metasources"]
    if 'metastatus' in status_group:
        metadata_group = status_group.metastatus
    else:
        metadata_group = h5f.create_group(status_group, "metastatus", "Statuses of the pangenome metadata")
    if 'metasources' in status_group:
        metasources_group = status_group.metasources
    else:
        metasources_group = h5f.create_group(status_group, "metasources", "Sources of the pangenome metadata")
    if metastatus["families"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.families = True
        metasources_group._v_attrs.families = metasources["families"]
    if metastatus["genes"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.genes = True
        metasources_group._v_attrs.genes = metasources["genes"]
    if metastatus["genomes"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.genomes = True
        metasources_group._v_attrs.genomes = metasources["genomes"]
    if metastatus["RGPs"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.RGPs = True
        metasources_group._v_attrs.RGPs = metasources["genes"]
    if metastatus["spots"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.spots = True
        metasources_group._v_attrs.spots = metasources["genes"]
    if metastatus["modules"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.modules = True
        metasources_group._v_attrs.modules = metasources["genes"]
    return True if any(metadata_group._v_attrs._f_list()) else False


def write_metadata_group(h5f: tables.File, metatype: str):
    if '/metadata' not in h5f:
        metadata_group = h5f.create_group("/", "metadata", "Pangenome metadata")
    else:
        metadata_group = h5f.root.metadata
    if metatype not in metadata_group:
        metatype_group = h5f.create_group(metadata_group, metatype, f"{metatype} metadata")
    else:
        metatype_group = metadata_group._f_get_child(metatype)
    return metatype_group


def write_row(meta_row, metadata):
    meta_row["value"] = metadata.value
    meta_row["accession"] = metadata.accession
    meta_row["description"] = metadata.description
    meta_row["score"] = metadata.score
    meta_row["e_val"] = metadata.e_val
    meta_row["bias"] = metadata.bias


def metadata_desc(max_name_len: int = 1, max_accession_len: int = 1, max_description_len: int = 1) -> dict:
    """
    Create a formated table for gene families metadata description
    :param max_name_len:
    :param max_accession_len:
    :param max_description_len:

    :return: Formated table
    """
    return {
        "accession": tables.StringCol(itemsize=max_accession_len),
        "value": tables.StringCol(itemsize=max_name_len),
        "description": tables.StringCol(itemsize=max_description_len),
        "score": tables.Float64Col(),
        "e_val": tables.Float64Col(),
        "bias": tables.Float64Col(),
    }


def desc_metadata_families(max_name_len: int = 1, max_accession_len: int = 1,
                           max_secondary_names_len: int = 1, max_description_len: int = 1,
                           max_gf_name_len: int = 1) -> dict:
    meta_desc = metadata_desc(max_name_len, max_accession_len, max_description_len)
    meta_desc.update({"secondary_names": tables.StringCol(itemsize=max_secondary_names_len),
                      "geneFam": tables.StringCol(itemsize=max_gf_name_len)})
    return meta_desc


def get_metadata_families_len(select_gf: List[GeneFamily], source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_gf: selected gene families from source
    :param source: Name of the metadata source
    :return: Maximum size of each element
    """
    max_value_len, max_accession_len, max_secondary_names_len, max_description_len = (1, 1, 1, 1)
    max_gf_name_len, expected_rows = (1, 0)

    for gf in select_gf:
        if len(gf.name) > max_gf_name_len:
            max_gf_name_len = len(gf.name)
        for metadata in gf.get_source(name=source):
            if len(metadata.value) > max_value_len:
                max_value_len = len(metadata.value)
            if metadata.accession is not None and len(metadata.accession) > max_accession_len:
                max_accession_len = len(metadata.accession)
            if metadata.description is not None and len(metadata.description) > max_description_len:
                max_description_len = len(metadata.description)
            if metadata.secondary_names is not None and len(metadata.secondary_names) > max_secondary_names_len:
                max_secondary_names_len = len(metadata.secondary_names)
            expected_rows += 1

    return max_value_len, max_accession_len, max_secondary_names_len, max_description_len, max_gf_name_len, expected_rows


def write_metadata_families(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, "families")
    select_gf = list(pangenome.get_gf_by_sources(source=source))
    meta_len = get_metadata_families_len(select_gf, source)
    source_table = h5f.create_table(metatype_group, source, desc_metadata_families(*meta_len[:-1]),
                                    expectedrows=meta_len[-1])
    meta_row = source_table.row
    for gf in tqdm(select_gf, unit='Gene family', desc=f'Source = {source}', disable=disable_bar):
        for metadata in gf.get_source(name=source):
            write_row(meta_row, metadata)
            meta_row["secondary_names"] = metadata.secondary_names
            meta_row["geneFam"] = gf.name
            meta_row.append()
    source_table.flush()


def desc_metadata_org(max_name_len: int = 1, max_accession_len: int = 1, max_description_len: int = 1,
                      max_org_name_len: int = 1) -> dict:
    meta_desc = metadata_desc(max_name_len, max_accession_len, max_description_len)
    meta_desc.update({"genome": tables.StringCol(itemsize=max_org_name_len)})
    return meta_desc


def get_metadata_org_len(select_org: List[Organism], source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_org: selected gene families from source
    :param source: Name of the metadata source
    :return: Maximum size of each element
    """
    max_value_len, max_accession_len, max_description_len = (1, 1, 1)
    max_org_name_len, expected_rows = (1, 0)

    for org in select_org:
        if len(org.name) > max_org_name_len:
            max_org_name_len = len(org.name)
        for metadata in org.get_source(name=source):
            if len(metadata.value) > max_value_len:
                max_value_len = len(metadata.value)
            if metadata.accession is not None and len(metadata.accession) > max_accession_len:
                max_accession_len = len(metadata.accession)
            if metadata.description is not None and len(metadata.description) > max_description_len:
                max_description_len = len(metadata.description)
            expected_rows += 1

    return max_value_len, max_accession_len, max_description_len, max_org_name_len, expected_rows


def write_metadata_genomes(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, "genomes")
    select_genomes = list(pangenome.get_org_by_sources(source=source))
    meta_len = get_metadata_org_len(select_genomes, source)
    source_table = h5f.create_table(metatype_group, source, desc_metadata_org(*meta_len[:-1]),
                                    expectedrows=meta_len[-1])
    meta_row = source_table.row
    for genome in tqdm(select_genomes, unit='Genome', desc=f'Source = {source}', disable=disable_bar):
        for metadata in genome.get_source(name=source):
            write_row(meta_row, metadata)
            meta_row["genome"] = genome.name
            meta_row.append()
    source_table.flush()


def desc_metadata_gene(max_name_len: int = 1, max_accession_len: int = 1, max_description_len: int = 1,
                       max_gene_name_len: int = 1) -> dict:
    meta_desc = metadata_desc(max_name_len, max_accession_len, max_description_len)
    meta_desc.update({"gene": tables.StringCol(itemsize=max_gene_name_len)})
    return meta_desc


def get_metadata_gene_len(select_gene: List[Gene], source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_gene: selected gene families from source
    :param source: Name of the metadata source
    :return: Maximum size of each element
    """
    max_value_len, max_accession_len, max_description_len = (1, 1, 1)
    max_gene_name_len, expected_rows = (1, 0)

    for gene in select_gene:
        if len(gene.name) > max_gene_name_len:
            max_gene_name_len = len(gene.name)
        for metadata in gene.get_source(name=source):
            if len(metadata.value) > max_value_len:
                max_value_len = len(metadata.value)
            if metadata.accession is not None and len(metadata.accession) > max_accession_len:
                max_accession_len = len(metadata.accession)
            if metadata.description is not None and len(metadata.description) > max_description_len:
                max_description_len = len(metadata.description)
            expected_rows += 1

    return max_value_len, max_accession_len, max_description_len, max_gene_name_len, expected_rows


def write_metadata_genes(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, "gene")
    select_gene = list(pangenome.get_gene_by_sources(source=source))
    meta_len = get_metadata_gene_len(select_gene, source)
    source_table = h5f.create_table(metatype_group, source, desc_metadata_gene(*meta_len[:-1]),
                                    expectedrows=meta_len[-1])
    meta_row = source_table.row
    for gene in tqdm(select_gene, unit='Gene', desc=f'Source = {source}', disable=disable_bar):
        for metadata in gene.get_source(name=source):
            write_row(meta_row, metadata)
            meta_row["gene"] = gene.name
            meta_row.append()
    source_table.flush()


def desc_metadata_rgp(max_name_len: int = 1, max_accession_len: int = 1, max_description_len: int = 1,
                      max_gene_name_len: int = 1) -> dict:
    meta_desc = metadata_desc(max_name_len, max_accession_len, max_description_len)
    meta_desc.update({"rgp": tables.StringCol(itemsize=max_gene_name_len)})
    return meta_desc


def get_metadata_rgp_len(select_rgps: List[Region], source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_rgps: selected gene families from source
    :param source: Name of the metadata source
    :return: Maximum size of each element
    """
    max_value_len, max_accession_len, max_description_len = (1, 1, 1)
    max_rgp_name_len, expected_rows = (1, 0)

    for rgp in select_rgps:
        if len(rgp.name) > max_rgp_name_len:
            max_rgp_name_len = len(rgp.name)
        for metadata in rgp.get_source(name=source):
            if len(metadata.value) > max_value_len:
                max_value_len = len(metadata.value)
            if metadata.accession is not None and len(metadata.accession) > max_accession_len:
                max_accession_len = len(metadata.accession)
            if metadata.description is not None and len(metadata.description) > max_description_len:
                max_description_len = len(metadata.description)
            expected_rows += 1

    return max_value_len, max_accession_len, max_description_len, max_rgp_name_len, expected_rows


def write_metadata_rgps(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, "RGP")
    select_rgps = list(pangenome.get_rgp_by_sources(source=source))
    meta_len = get_metadata_rgp_len(select_rgps, source)
    source_table = h5f.create_table(metatype_group, source, desc_metadata_rgp(*meta_len[:-1]),
                                    expectedrows=meta_len[-1])
    meta_row = source_table.row
    for rgp in tqdm(select_rgps, unit='RGP', desc=f'Source = {source}', disable=disable_bar):
        for metadata in rgp.get_source(name=source):
            write_row(meta_row, metadata)
            meta_row["rgp"] = rgp.name
            meta_row.append()
    source_table.flush()


def desc_metadata_spot(max_name_len: int = 1, max_accession_len: int = 1, max_description_len: int = 1,
                       max_spot_name_len: int = 1) -> dict:
    meta_desc = metadata_desc(max_name_len, max_accession_len, max_description_len)
    meta_desc.update({"spot": tables.StringCol(itemsize=max_spot_name_len)})
    return meta_desc


def get_metadata_spots_len(select_spots: List[Spot], source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_spots: selected gene families from source
    :param source: Name of the metadata source
    :return: Maximum size of each element
    """
    max_value_len, max_accession_len, max_description_len = (1, 1, 1)
    max_spot_name_len, expected_rows = (1, 0)

    for spot in select_spots:
        if len(spot.ID) > max_spot_name_len:
            max_spot_name_len = len(spot.ID)
        for metadata in spot.get_source(name=source):
            if len(metadata.value) > max_value_len:
                max_value_len = len(metadata.value)
            if metadata.accession is not None and len(metadata.accession) > max_accession_len:
                max_accession_len = len(metadata.accession)
            if metadata.description is not None and len(metadata.description) > max_description_len:
                max_description_len = len(metadata.description)
            expected_rows += 1

    return max_value_len, max_accession_len, max_description_len, max_spot_name_len, expected_rows


def write_metadata_spots(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, "spot")
    select_spots = list(pangenome.get_spots_by_sources(source=source))
    meta_len = get_metadata_spots_len(select_spots, source)
    source_table = h5f.create_table(metatype_group, source, desc_metadata_spot(*meta_len[:-1]),
                                    expectedrows=meta_len[-1])
    meta_row = source_table.row
    for spot in tqdm(select_spots, unit='Spot', desc=f'Source = {source}', disable=disable_bar):
        for metadata in spot.get_source(name=source):
            write_row(meta_row, metadata)
            meta_row["spot"] = spot.ID
            meta_row.append()
    source_table.flush()


def desc_metadata_modules(max_name_len: int = 1, max_accession_len: int = 1, max_description_len: int = 1,
                         max_module_name_len: int = 1) -> dict:
    meta_desc = metadata_desc(max_name_len, max_accession_len, max_description_len)
    meta_desc.update({"module": tables.StringCol(itemsize=max_module_name_len)})
    return meta_desc


def get_metadata_modules_len(select_modules: List[Module], source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_modules: selected gene families from source
    :param source: Name of the metadata source
    :return: Maximum size of each element
    """
    max_value_len, max_accession_len, max_description_len = (1, 1, 1)
    max_module_name_len, expected_rows = (1, 0)

    for module in select_modules:
        if len(module.ID) > max_module_name_len:
            max_module_name_len = len(module.ID)
        for metadata in module.get_source(name=source):
            if len(metadata.value) > max_value_len:
                max_value_len = len(metadata.value)
            if metadata.accession is not None and len(metadata.accession) > max_accession_len:
                max_accession_len = len(metadata.accession)
            if metadata.description is not None and len(metadata.description) > max_description_len:
                max_description_len = len(metadata.description)
            expected_rows += 1

    return max_value_len, max_accession_len, max_description_len, max_module_name_len, expected_rows


def write_metadata_modules(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, "module")
    select_modules = list(pangenome.get_spots_by_sources(source=source))
    meta_len = get_metadata_modules_len(select_modules, source)
    source_table = h5f.create_table(metatype_group, source, desc_metadata_modules(*meta_len[:-1]),
                                    expectedrows=meta_len[-1])
    meta_row = source_table.row
    for module in tqdm(select_modules, unit='Module', desc=f'Source = {source}', disable=disable_bar):
        for metadata in module.get_source(name=source):
            write_row(meta_row, metadata)
            meta_row["module"] = module.ID
            meta_row.append()
    source_table.flush()


def erase_metadata(pangenome: Pangenome, h5f: tables.File, status_group: tables.Group,
                   metatype: str = None, source: str = None):
    metadata_group = h5f.root.metadata

    metastatus = pangenome.status["metadata"]
    metasources = pangenome.status["metasources"]
    if metatype in metadata_group:
        metatype_group = metadata_group._f_get_child(metatype)
        if source in metatype_group:
            logging.getLogger().info(f"Erasing metadata assign to {metatype} from source {source}...")
            metasources[metatype].remove(source)
            status_group.metasources._v_attrs[metatype].remove(source)
            h5f.remove_node(metatype_group, source)
        if metatype_group._v_nchildren == 0:
            logging.getLogger().debug(f"No more source of metadata in {metatype}."
                                      f"Erasing node {metatype} in metadata")
            metastatus[metatype] = 'No'
            status_group.metastatus.families = False
            h5f.remove_node(metadata_group, metatype)
    if metadata_group._v_nchildren == 0:
        logging.getLogger().debug("No more metadata in pangenome. Erasing node metadata")
        status_group._v_attrs.metadata = False
        h5f.remove_node("/", "metadata")
        h5f.remove_node(status_group, "metasources")
        h5f.remove_node(status_group, "metastatus")


def write_metadata(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    if pangenome.status["metadata"]["families"] == "Computed":
        logging.getLogger().info("Writing gene families metadata in pangenome")
        write_metadata_families(pangenome, h5f, pangenome.status["metasources"]["families"][-1], disable_bar)
        pangenome.status["metadata"]["families"] = "Loaded"

    if pangenome.status["metadata"]["genomes"] == "Computed":
        logging.getLogger().info("Writing genomes metadata in pangenome")
        write_metadata_genomes(pangenome, h5f, pangenome.status["metasources"]["genomes"][-1], disable_bar)
        pangenome.status["metadata"]["genomes"] = "Loaded"

    if pangenome.status["metadata"]["genes"] == "Computed":
        logging.getLogger().info("Writing genes metadata in pangenome")
        write_metadata_genes(pangenome, h5f, pangenome.status["metasources"]["genes"][-1], disable_bar)
        pangenome.status["metadata"]["genes"] = "Loaded"

    if pangenome.status["metadata"]["rgps"] == "Computed":
        logging.getLogger().info("Writing genes metadata in pangenome")
        write_metadata_rgps(pangenome, h5f, pangenome.status["metasources"]["rgps"][-1], disable_bar)
        pangenome.status["metadata"]["rgps"] = "Loaded"

    if pangenome.status["metadata"]["spots"] == "Computed":
        logging.getLogger().info("Writing genes metadata in pangenome")
        write_metadata_spots(pangenome, h5f, pangenome.status["metasources"]["spots"][-1], disable_bar)
        pangenome.status["metadata"]["spots"] = "Loaded"

    if pangenome.status["metadata"]["modules"] == "Computed":
        logging.getLogger().info("Writing genes metadata in pangenome")
        write_metadata_genomes(pangenome, h5f, pangenome.status["metasources"]["modules"][-1], disable_bar)
        pangenome.status["metadata"]["modules"] = "Loaded"
