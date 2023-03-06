#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import tempfile
import subprocess
import argparse
from collections import defaultdict
import sys
import pdb

# local libraries
from ppanggolin.formats import check_pangenome_info, write_pangenome
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import write_pangenome


def read_tsv_file_metadata(tsv_path):
    with open(tsv_path, 'r') as tsv_file:
        metadata_names = []
        ret_metadata = defaultdict(dict)
        for num, line in enumerate(tsv_file):
            elements = [e.strip() for e in line.split("\t")]
            if num == 0:
                metadata_names = elements[1:]
            else:
                for index_metadata_name, metadata_value in enumerate(elements[1:]):
                    ret_metadata[elements[0]][metadata_names[index_metadata_name]] = metadata_value
    return (ret_metadata)


def launch(args):
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    check_pangenome_info(pangenome, need_families=True, need_annotations=True, need_metadata=True)
    if args.on_families is not None:
        for fam_name, metadata_fam in read_tsv_file_metadata(args.on_families).items():
            pangenome.get_gene_family(fam_name).metadata.update(metadata_fam)
        pangenome.status["metadata_on_families"] = "Computed"
    if args.on_organisms is not None:
        for org_name, metadata_org in read_tsv_file_metadata(args.on_organisms).items():
            pangenome.get_organism(org_name).metadata.update(metadata_org)
        pangenome.status["metadata_on_organims"] = "Computed"
    write_pangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)


def metadataSubparser(subparser):
    parser = subparser.add_parser("metadata", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Mandatory input/output files",
                                         description="One of the following argument is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")

    parser.add_argument('--on_families', default=None, required=False, type=str,
                        help="import a tsv file containing metadata on families")
    parser.add_argument('--on_organisms', default=None, required=False, type=str,
                        help="import a tsv file containing metadata on organisms")

    return parser