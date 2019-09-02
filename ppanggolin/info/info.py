#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse

#installed libraries
import tables

#local libraries
from ppanggolin.formats import readInfo

def printInfo(args):
    h5f = tables.open_file(args.pangenome,"r")
    statusGroup = h5f.root.status
    print(f"genomes annotated : {'true' if statusGroup._v_attrs.genomesAnnotated else 'false' }")
    print(f"genes clustered : {'true' if statusGroup._v_attrs.genesClustered else 'false' }")
    print(f"genes have their sequences : {'true' if statusGroup._v_attrs.geneSequences else 'false' }")
    print(f"gene families have their sequences : {'true' if statusGroup._v_attrs.geneFamilySequences else 'false' }")
    print(f"neighbors graph : {'true' if statusGroup._v_attrs.NeighborsGraph else 'false' }")
    print(f"pangenome partitioned : {'true' if statusGroup._v_attrs.Partitionned else 'false' }")

    if args.details:
        readInfo(h5f)

    h5f.close()


def infoSubparser(subparser):
    parser = subparser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    options = parser.add_argument_group(title = "optional arguments")
    options.add_argument("--details",required=False,action="store_true", help="Shows detailled informations about the pangenome's content")
    required = parser.add_argument_group(title = "Required arguments", description = "The following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser
