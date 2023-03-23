#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from pathlib import Path
import csv

# installed libraries
from tqdm import tqdm
import pandas as pd

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.metadata import Metadata
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome

gf_meta_col_names = ['Gene_family', 'Accession', 'protein_name', 'e_value',
                     'score', 'bias', 'secondary_name', 'Description']
gf_meta_col_type = {'Gene_family': str,
                    'Accession': str,
                    'protein_name': str,
                    'e_value': float,
                    'score': float,
                    'bias': float,
                    'secondary_name': str,
                    'Description': str}


def check_pangenome_metadata(pangenome: Pangenome, source: str, force: bool = False, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param source: source of the annotation
    :param force: erase if an annotation for the provide source already exist
    :param disable_bar: Disable bar
    """
    if pangenome.status["metadata"]["families"] == "inFile" and source in pangenome.status["metasources"]["families"]:
        if force:
            erase_pangenome(pangenome, metadata=True, source=source, metatype="families")
        else:
            raise Exception(f"An annotation corresponding to the source : '{source}' already exist in pangenome."
                            "Add the option --force to erase")
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=disable_bar)


def assign_metadata_to_families(annotation_df: pd.DataFrame, pangenome: Pangenome, source: str = None,
                                omit: bool = False, disable_bar: bool = False) -> dict:
    """ Add to gene families an annotation and create a dictionary with for each annotation a set of gene family

    :param annotation_df: Dataframe with for each family an annotation
    :param pangenome: Pangenome with gene families
    :param source: source of the annotation
    :param disable_bar:
    :return: Dictionary with for each annotation a set of gene family
    """
    for row in tqdm(annotation_df.itertuples(index=False), unit='row',
                    total=annotation_df.shape[0], disable=disable_bar):
        try:
            gene_fam = pangenome.get_gene_family(name=row.Gene_family)
        except KeyError:
            if omit:
                pass
            else:
                raise KeyError(f"Family {row.Gene_family} does not exist in pangenome. Check name in your file")
        else:
            annotation = Metadata(source=source, accession=row.Accession, value=row.protein_name,
                                  secondary_names=row.secondary_name, description=row.Description,
                                  score=row.score, e_val=row.e_value, bias=row.bias)
            gene_fam.add_metadata(source=source, metadata=annotation)

    pangenome.status["metadata"]["families"] = "Computed"
    pangenome.status["metasources"]["families"].append(source)


def metadata_to_families(pangenome: Pangenome, data: Path, source: str, omit: bool = False,
                         force: bool = False, disable_bar: bool = False):
    try:
        assert "metadata" in pangenome.status and "metasources" in pangenome.status
        assert "families" in pangenome.status["metadata"] and "families" in pangenome.status["metasources"]
    except AssertionError:
        raise AssertionError("Unexpected problem. Please report the error in our GitHub")

    check_pangenome_metadata(pangenome, source=source, force=force, disable_bar=disable_bar)
    annotation_df = pd.read_csv(data, sep="\t", header=None, quoting=csv.QUOTE_NONE,
                                names=gf_meta_col_names, dtype=gf_meta_col_type)
    assign_metadata_to_families(annotation_df, pangenome, source, omit, disable_bar)
    logging.getLogger().info("Metadata assignment Done")
    write_pangenome(pangenome, pangenome.file, disable_bar=disable_bar)
