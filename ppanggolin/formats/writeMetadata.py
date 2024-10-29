#!/usr/bin/env python3

# default libraries
import logging
from typing import Dict, List, Tuple, Union

# installed libraries
from tqdm import tqdm
import numpy
import tables

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Organism, Gene, Contig
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot, Module


def write_metadata_status(
    pangenome: Pangenome, h5f: tables.File, status_group: tables.Group
) -> bool:
    """Write status of metadata in pangenome file

    :param pangenome: pangenome with metadata
    :param h5f: HDF5 file with pangenome
    :param status_group: Pangenome status information group
    :return:
    """
    metastatus = pangenome.status["metadata"]
    metasources = pangenome.status["metasources"]
    if "metastatus" in status_group:
        metadata_group = status_group.metastatus
    else:
        metadata_group = h5f.create_group(
            status_group, "metastatus", "Statuses of the pangenome metadata"
        )
    if "metasources" in status_group:
        metasources_group = status_group.metasources
    else:
        metasources_group = h5f.create_group(
            status_group, "metasources", "Sources of the pangenome metadata"
        )

    if metastatus["families"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.families = True
        metasources_group._v_attrs.families = metasources["families"]
    if metastatus["genes"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.genes = True
        metasources_group._v_attrs.genes = metasources["genes"]
    if metastatus["contigs"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.contigs = True
        metasources_group._v_attrs.contigs = metasources["contigs"]
    if metastatus["genomes"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.genomes = True
        metasources_group._v_attrs.genomes = metasources["genomes"]
    if metastatus["RGPs"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.RGPs = True
        metasources_group._v_attrs.RGPs = metasources["RGPs"]
    if metastatus["spots"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.spots = True
        metasources_group._v_attrs.spots = metasources["spots"]
    if metastatus["modules"] in ["Computed", "Loaded", "inFile"]:
        metadata_group._v_attrs.modules = True
        metasources_group._v_attrs.modules = metasources["modules"]

    return any(metadata_group._v_attrs._f_list())


def write_metadata_group(h5f: tables.File, metatype: str) -> tables.Group:
    """Check and write the group in HDF5 file to organize metadata

    :param h5f: HDF5 file with pangenome
    :param metatype: select to which pangenome element metadata should be written

    :return: Metadata group of the corresponding metatype
    """
    if "/metadata" not in h5f:
        metadata_group = h5f.create_group("/", "metadata", "Pangenome metadata")
    else:
        metadata_group = h5f.root.metadata
    if metatype not in metadata_group:
        metatype_group = h5f.create_group(
            metadata_group, metatype, f"{metatype} metadata"
        )
    else:
        metatype_group = metadata_group._f_get_child(metatype)
    return metatype_group


def desc_metadata(
    max_len_dict: Dict[str, int], type_dict: Dict[str, tables.Col]
) -> dict:
    """Create a formatted table for metadata description

    :return: Formatted table
    """
    desc_dict = {
        attr: tables.StringCol(itemsize=max_value)
        for attr, max_value in max_len_dict.items()
    }
    desc_dict.update(dict(type_dict.items()))
    return desc_dict


def get_metadata_contig_len(
    select_ctg: List[Contig], source: str
) -> Tuple[Dict[str, int], Dict[str, tables.Col], int]:
    """Get maximum size of contig metadata information

    :param select_ctg: selected elements from source
    :param source: Name of the metadata source

    :return: Maximum type and size of each element
    """
    type_dict = {"metadata_id": tables.Int64Col(), "ID": tables.Int64Col()}
    max_len_dict = {}
    expected_rows = 0

    for contig in select_ctg:
        for metadata in contig.get_metadata_by_source(source).values():
            for attr, value in (
                (k, v) for k, v in metadata.__dict__.items() if k != "source"
            ):
                if isinstance(value, bytes):
                    value = value.decode("UTF-8")
                if isinstance(value, float) or isinstance(value, int):
                    if attr in type_dict:
                        if type_dict[attr] != type(value):
                            if isinstance(value, float) and type_dict[attr] == int:
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
                    logging.getLogger("PPanGGOLiN").debug(
                        f"attr: {attr}, value: {value}"
                    )
                    raise TypeError(f"{type(value)} is not acceptable")
            expected_rows += 1

    return max_len_dict, type_dict, expected_rows


def write_metadata_contig(
    h5f: tables.File,
    source: str,
    select_contigs: List[Contig],
    disable_bar: bool = False,
):
    """Writing a table containing the metadata associated to contig

    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param select_contigs: List of contig withj metadata
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, "contigs")
    meta_len = get_metadata_contig_len(select_contigs, source)
    # h5f.remove_node(metatype_group, source)
    source_table = h5f.create_table(
        metatype_group, source, desc_metadata(*meta_len[:-1]), expectedrows=meta_len[-1]
    )
    meta_row = source_table.row
    for contig in tqdm(
        select_contigs, unit="contigs", desc=f"Source = {source}", disable=disable_bar
    ):
        for meta_id, metadata in contig.get_metadata_by_source(source).items():
            meta_row["metadata_id"] = meta_id
            for desc in source_table.colnames:
                if desc == "ID":
                    meta_row[desc] = contig.ID
                elif desc == "metadata_id":
                    meta_row[desc] = meta_id
                else:
                    if hasattr(metadata, desc):
                        value = metadata.__getattribute__(desc)
                        if isinstance(value, bytes):
                            value = value.decode("UTF-8")
                        meta_row[desc] = value
                    else:
                        meta_row[desc] = None
            meta_row.append()
    source_table.flush()


def get_metadata_len(
    select_elem: Union[
        List[Gene],
        List[Organism],
        List[GeneFamily],
        List[Region],
        List[Spot],
        List[Module],
    ],
    source: str,
) -> Tuple[Dict[str, int], Dict[str, tables.Col], int]:
    """Get maximum size of metadata information

    :param select_elem: selected elements from source
    :param source: Name of the metadata source

    :return: Maximum type and size of each element
    """
    type_dict = {"metadata_id": tables.Int64Col()}
    max_len_dict = {}
    expected_rows = 0

    for element in select_elem:
        if hasattr(element, "name") and element.name:
            max_len_dict["ID"] = max(max_len_dict.get("ID", 0), len(element.name))
        elif hasattr(element, "ID"):
            if isinstance(element.ID, str):
                max_len_dict["ID"] = max(max_len_dict.get("ID", 0), len(element.ID))
            elif isinstance(
                element.ID, (int, numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64)
            ):
                type_dict["ID"] = tables.Int64Col()
            else:
                raise TypeError(
                    f"Invalid type for 'ID' in element '{element}': expected integer-like type but got "
                    f"{type(element.ID).__name__}."
                )
        else:
            raise AttributeError(
                f"Unexpected attribute in element '{element}': missing 'name' or 'ID'. "
                "Please report this error on our GitHub."
            )

        for metadata in element.get_metadata_by_source(source).values():
            for attr, value in (
                (k, v) for k, v in metadata.__dict__.items() if k != "source"
            ):
                if isinstance(value, bytes):
                    value = value.decode("UTF-8")
                if isinstance(value, float) or isinstance(value, int):
                    if attr in type_dict:
                        if isinstance(type_dict[attr], type(value)):
                            if isinstance(value, float) and isinstance(
                                type_dict[attr], int
                            ):
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
                    logging.getLogger("PPanGGOLiN").debug(
                        f"attr: {attr}, value: {value}"
                    )
                    raise TypeError(
                        f"Invalid metadata type: The attribute '{attr}' from the pangenome element '{element}' "
                        f"has an unexpected value '{value}' of type '{type(value).__name__}'."
                    )

            expected_rows += 1

    for attribute, max_length in max_len_dict.items():
        if max_length == 0:
            raise ValueError(
                f"Metadata attribute '{attribute}' has a length of 0, which is not allowed."
            )

    return max_len_dict, type_dict, expected_rows


def write_metadata_metatype(
    h5f: tables.File,
    source: str,
    metatype: str,
    select_elements: Union[
        List[Gene],
        List[Organism],
        List[GeneFamily],
        List[Region],
        List[Spot],
        List[Module],
    ],
    disable_bar: bool = False,
):
    """Writing a table containing the metadata associated to element from the metatype

    :param h5f: HDF5 file to write gene families
    :param source: name of the metadata source
    :param metatype: select to which pangenome element metadata should be written
    :param select_elements: Elements selected to write metadata
    :param disable_bar: Disable progress bar
    """
    metatype_group = write_metadata_group(h5f, metatype)
    max_len_dict, type_dict, expected_rows = get_metadata_len(select_elements, source)

    desc_metadata(max_len_dict, type_dict)

    source_table = h5f.create_table(
        metatype_group,
        source,
        desc_metadata(max_len_dict, type_dict),
        expectedrows=expected_rows,
    )
    meta_row = source_table.row
    for element in tqdm(
        select_elements, unit=metatype, desc=f"Source = {source}", disable=disable_bar
    ):
        for meta_id, metadata in element.get_metadata_by_source(source).items():
            for desc in source_table.colnames:
                if desc == "ID":
                    if hasattr(element, "name") and len(element.name) > 0:
                        meta_row[desc] = element.name
                    else:
                        meta_row[desc] = element.ID
                elif desc == "metadata_id":
                    meta_row[desc] = meta_id
                else:
                    if hasattr(metadata, desc):
                        value = metadata.__getattribute__(desc)
                        if isinstance(value, bytes):
                            value = value.decode("UTF-8")
                        meta_row[desc] = value
            meta_row.append()
    source_table.flush()


def erase_metadata(
    pangenome: Pangenome,
    h5f: tables.File,
    status_group: tables.Group,
    metatype: str = None,
    source: str = None,
):
    """
    Erase metadata in pangenome

    :param pangenome: Pangenome with metadata to erase
    :param h5f: HDF5 file with pangenome metadata
    :param status_group: pangenome status in HDF5
    :param metatype: select to which pangenome element metadata should be erased
    :param source: name of the metadata source
    """
    metadata_group = h5f.root.metadata

    metastatus = pangenome.status["metadata"]
    metasources = pangenome.status["metasources"]
    if metatype in metadata_group:
        metatype_group = metadata_group._f_get_child(metatype)
        if source in metatype_group:
            logging.getLogger("PPanGGOLiN").info(
                f"Erasing metadata assign to {metatype} from source {source}..."
            )
            metasources[metatype].remove(source)
            status_group.metasources._v_attrs[metatype].remove(source)
            h5f.remove_node(metatype_group, source)
        if metatype_group._v_nchildren == 0:
            logging.getLogger("PPanGGOLiN").debug(
                f"No more source of metadata in {metatype}. "
                f"Erasing node {metatype} in metadata"
            )
            metastatus[metatype] = "No"
            status_group.metastatus.families = False
            h5f.remove_node(metadata_group, metatype)
    if metadata_group._v_nchildren == 0:
        logging.getLogger("PPanGGOLiN").debug(
            "No more metadata in pangenome. Erasing node metadata"
        )
        status_group._v_attrs.metadata = False
        h5f.remove_node("/", "metadata")
        h5f.remove_node(status_group, "metasources")
        h5f.remove_node(status_group, "metastatus")


def write_metadata(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """Write metadata in pangenome

    :param pangenome: Pangenome where should be written metadata
    :param h5f: HDF5 file with pangenome
    :param disable_bar: Disable progress bar
    """
    if pangenome.status["metadata"]["families"] == "Computed":
        logging.getLogger("PPanGGOLiN").info(
            "Writing gene families metadata in pangenome"
        )
        select_gf = list(
            pangenome.get_elem_by_source(
                source=pangenome.status["metasources"]["families"][-1],
                metatype="families",
            )
        )
        write_metadata_metatype(
            h5f,
            pangenome.status["metasources"]["families"][-1],
            "families",
            select_gf,
            disable_bar,
        )
        pangenome.status["metadata"]["families"] = "Loaded"

    if pangenome.status["metadata"]["genomes"] == "Computed":
        logging.getLogger("PPanGGOLiN").info("Writing genomes metadata in pangenome")
        select_genomes = list(
            pangenome.get_elem_by_source(
                source=pangenome.status["metasources"]["genomes"][-1],
                metatype="genomes",
            )
        )
        write_metadata_metatype(
            h5f,
            pangenome.status["metasources"]["genomes"][-1],
            "genomes",
            select_genomes,
            disable_bar,
        )
        pangenome.status["metadata"]["genomes"] = "Loaded"

    if pangenome.status["metadata"]["contigs"] == "Computed":
        logging.getLogger("PPanGGOLiN").info("Writing contigs metadata in pangenome")
        select_contigs = list(
            pangenome.get_elem_by_source(
                source=pangenome.status["metasources"]["contigs"][-1],
                metatype="contigs",
            )
        )
        write_metadata_contig(
            h5f,
            pangenome.status["metasources"]["contigs"][-1],
            select_contigs,
            disable_bar,
        )
        pangenome.status["metadata"]["contigs"] = "Loaded"

    if pangenome.status["metadata"]["genes"] == "Computed":
        logging.getLogger("PPanGGOLiN").info("Writing genes metadata in pangenome")
        select_genes = list(
            pangenome.get_elem_by_source(
                source=pangenome.status["metasources"]["genes"][-1], metatype="genes"
            )
        )
        write_metadata_metatype(
            h5f,
            pangenome.status["metasources"]["genes"][-1],
            "genes",
            select_genes,
            disable_bar,
        )
        pangenome.status["metadata"]["genes"] = "Loaded"

    if pangenome.status["metadata"]["RGPs"] == "Computed":
        logging.getLogger("PPanGGOLiN").info("Writing RGPs metadata in pangenome")
        select_rgps = list(
            pangenome.get_elem_by_source(
                source=pangenome.status["metasources"]["RGPs"][-1], metatype="RGPs"
            )
        )
        write_metadata_metatype(
            h5f,
            pangenome.status["metasources"]["RGPs"][-1],
            "RGPs",
            select_rgps,
            disable_bar,
        )
        pangenome.status["metadata"]["RGPs"] = "Loaded"

    if pangenome.status["metadata"]["spots"] == "Computed":
        logging.getLogger("PPanGGOLiN").info("Writing spots metadata in pangenome")
        select_spots = list(
            pangenome.get_elem_by_source(
                source=pangenome.status["metasources"]["spots"][-1], metatype="spots"
            )
        )
        write_metadata_metatype(
            h5f,
            pangenome.status["metasources"]["spots"][-1],
            "spots",
            select_spots,
            disable_bar,
        )
        pangenome.status["metadata"]["spots"] = "Loaded"

    if pangenome.status["metadata"]["modules"] == "Computed":
        logging.getLogger("PPanGGOLiN").info("Writing modules metadata in pangenome")
        select_modules = list(
            pangenome.get_elem_by_source(
                source=pangenome.status["metasources"]["modules"][-1],
                metatype="modules",
            )
        )
        write_metadata_metatype(
            h5f,
            pangenome.status["metasources"]["modules"][-1],
            "modules",
            select_modules,
            disable_bar,
        )
        pangenome.status["metadata"]["modules"] = "Loaded"
