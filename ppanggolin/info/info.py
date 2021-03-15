#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse

#installed libraries
import tables

#local libraries
from ppanggolin.formats import readInfo, readParameters

def printInfo(pangenome, status = False, content = False, parameters = False):
    if status or content or parameters:
        h5f = tables.open_file(pangenome,"r")
        if status:
            statusGroup = h5f.root.status
            print(f"genomes annotated : {'true' if statusGroup._v_attrs.genomesAnnotated else 'false' }")
            print(f"genes clustered : {'true' if statusGroup._v_attrs.genesClustered else 'false' }")
            print(f"genes have their sequences : {'true' if statusGroup._v_attrs.geneSequences else 'false' }")
            print(f"gene families have their sequences : {'true' if statusGroup._v_attrs.geneFamilySequences else 'false' }")
            print(f"neighbors graph : {'true' if statusGroup._v_attrs.NeighborsGraph else 'false' }")
            print(f"pangenome partitioned : {'true' if statusGroup._v_attrs.Partitionned else 'false' }")
            if hasattr(statusGroup._v_attrs, "predictedRGP"):
                print(f"RGP predicted : {'true' if statusGroup._v_attrs.predictedRGP else 'false' }")

            if hasattr(statusGroup._v_attrs, "spots"):
                print(f"Spots predicted : {'true' if statusGroup._v_attrs.spots else 'false' }")

            if hasattr(statusGroup._v_attrs, "version"):
                print(f"PPanGGOLiN version : {statusGroup._v_attrs.version}")

        if content:
            readInfo(h5f)
        if parameters:
            readParameters(h5f)
        h5f.close()
    else:
        print("Please select what information you want by using --parameters, --content or --status")

def launch(args):
    printInfo(args.pangenome, args.status, args.content, args.parameters)

def infoSubparser(subparser):
    parser = subparser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Required arguments", description = "The following arguments is required :")
    required.add_argument('-p','--pangenome', required=True, type=str, help="The pangenome .h5 file")

    options = parser.add_argument_group(title = "optional arguments")
    options.add_argument("--parameters", required=False, action = "store_true", help = "Shows the parameters used (or computed) for each step of the pangenome generation")
    options.add_argument("--content", required=False, action="store_true", help="Shows detailled informations about the pangenome's content")
    options.add_argument("--status", required=False, action="store_true", help="Shows informations about the statuses of the different elements of the pangenome (what has been computed, or not)")

    return parser
