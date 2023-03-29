#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from typing import Dict, List, Tuple, Union

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


def desc_metadata(max_len_dict: Dict[str, int], type_dict: Dict[str, tables.Col]) -> dict:
    """
    Create a formated table for gene families metadata description

    :return: Formated table
    """
    desc_dict = {attr: tables.StringCol(itemsize=max_value) for attr, max_value in max_len_dict.items()}
    desc_dict.update({attr: col_type for attr, col_type in type_dict.items()})
    return desc_dict


def get_metadata_len(select_elem: List[Module], source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_elem: selected gene families from source
    :param source: Name of the metadata source
    :return: Maximum size of each element
    """
    type_dict = {}
    max_len_dict = {}
    expected_rows = 0

    for element in select_elem:
        if hasattr(element, 'name'):
            if 'ID' in max_len_dict:
                if len(element.name) > max_len_dict['ID']:
                    max_len_dict['ID'] = len(element.name)
            else:
                if "name" not in max_len_dict or len(element.name) > max_len_dict['name']:
                    max_len_dict['name'] = len(element.name)
        elif hasattr(element, 'ID'):
            if 'name' in max_len_dict:
                if len(element.ID) > max_len_dict['name']:
                    max_len_dict['name'] = len(element.ID)
        else:
            raise Exception("Unexpected attribute. A recent change could create this error."
                            " Please report the error on our github.")
        for metadata in element.get_source(name=source):
            for attr, value in ((k, v) for k, v in metadata.__dict__.items() if k != "source"):
                if isinstance(value, bytes):
                    value = value.decode('UTF-8')
                if isinstance(value, float) or isinstance(value, int):
                    if attr in type_dict:
                        if type_dict[attr] != type(value):
                            if type(value) == float and type_dict[attr] == int:
                                type_dict[attr] = tables.Float64Col()
                    else:
                        if isinstance(value, float):
                            type_dict[attr] = tables.Float64Col()
                        else:
                            type_dict[attr] = tables.Int64Col()
                elif isinstance(value, str):
                    if attr in max_len_dict:
                        if len(value) > max_len_dict[attr]:
                            max_len_dict[attr] = len(value)
                    else:
                        max_len_dict[attr] = len(value)
                else:
                    raise TypeError(f"{type(value)} is not acceptable")
            expected_rows += 1

    return max_len_dict, type_dict, expected_rows


def write_metadata_metatype(h5f: tables.File, source: str, metatype: str,
                            select_elements: Union[List[Gene], List[Organism], List[GeneFamily], List[Region],
                                             List[Spot], List[Module]],
                            disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, metatype)
    meta_len = get_metadata_len(select_elements, source)
    source_table = h5f.create_table(metatype_group, source, desc_metadata(*meta_len[:-1]), expectedrows=meta_len[-1])
    meta_row = source_table.row
    for element in tqdm(select_elements, unit=metatype, desc=f'Source = {source}', disable=disable_bar):
        for metadata in element.get_source(name=source):
            for desc in source_table.colnames:
                if desc in ["name", "ID"]:
                    meta_row[desc] = element.__getattribute__(desc)
                else:
                    value = metadata.__getattribute__(desc)
                    if isinstance(value, bytes):
                        value = value.decode('UTF-8')
                    meta_row[desc] = value
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
        select_gf = list(pangenome.get_gf_by_sources(source=pangenome.status["metasources"]["families"][-1]))
        write_metadata_metatype(h5f, pangenome.status["metasources"]["families"][-1],
                                "families", select_gf, disable_bar)
        pangenome.status["metadata"]["families"] = "Loaded"

    if pangenome.status["metadata"]["genomes"] == "Computed":
        logging.getLogger().info("Writing genomes metadata in pangenome")
        select_genomes = list(pangenome.get_org_by_sources(source=pangenome.status["metasources"]["genomes"][-1]))
        write_metadata_metatype(h5f, pangenome.status["metasources"]["genomes"][-1],
                                "genomes", select_genomes, disable_bar)
        pangenome.status["metadata"]["genomes"] = "Loaded"

    if pangenome.status["metadata"]["genes"] == "Computed":
        logging.getLogger().info("Writing genes metadata in pangenome")
        select_genes = list(pangenome.get_gene_by_sources(source=pangenome.status["metasources"]["genes"][-1]))
        write_metadata_metatype(h5f, pangenome.status["metasources"]["genes"][-1],
                                "genes", select_genes, disable_bar)
        pangenome.status["metadata"]["genes"] = "Loaded"

    if pangenome.status["metadata"]["RGPs"] == "Computed":
        logging.getLogger().info("Writing genes metadata in pangenome")
        select_rgps = list(pangenome.get_rgp_by_sources(source=pangenome.status["metasources"]["RGPs"][-1]))
        write_metadata_metatype(h5f, pangenome.status["metasources"]["RGPs"][-1],
                                "RGPs", select_rgps, disable_bar)
        pangenome.status["metadata"]["RGPs"] = "Loaded"

    if pangenome.status["metadata"]["spots"] == "Computed":
        logging.getLogger().info("Writing genes metadata in pangenome")
        select_spots = list(pangenome.get_spots_by_sources(source=pangenome.status["metasources"]["spots"][-1]))
        write_metadata_metatype(h5f, pangenome.status["metasources"]["spots"][-1],
                                "spots", select_spots, disable_bar)
        pangenome.status["metadata"]["spots"] = "Loaded"

    if pangenome.status["metadata"]["modules"] == "Computed":
        logging.getLogger().info("Writing genes metadata in pangenome")
        select_modules = list(pangenome.get_modules_by_sources(source=pangenome.status["metasources"]["modules"][-1]))
        write_metadata_metatype(h5f, pangenome.status["metasources"]["modules"][-1],
                                "modules", select_modules, disable_bar)
        pangenome.status["metadata"]["modules"] = "Loaded"
