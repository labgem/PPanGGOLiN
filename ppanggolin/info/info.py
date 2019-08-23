#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse

#installed libraries
import tables

def printInfo(pangenomeFile):
    h5f = tables.open_file(pangenomeFile,"r")
    statusGroup = h5f.root.status
    print("genomes annotated : " + "true" if statusGroup._v_attrs.genomesAnnotated else "false")
    print("genes clustered : " + "true" if statusGroup._v_attrs.genesClustered else "false")
    print("genes have their sequences : " + "true" if statusGroup._v_attrs.geneSequences else "false")
    print("gene families have their sequences : " + "true" if statusGroup._v_attrs.geneFamilySequences else "false")
    print("neighbors graph : " + "true" if statusGroup._v_attrs.NeighborsGraph else "false")
    print("pangenome partitioned : " + "true" if statusGroup._v_attrs.Partitionned else "false")
    h5f.close()


def infoSubparser(subparser):
    parser = subparser.add_parser("info",help = "Prints information about a given pangenome graph file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title = "Required arguments", description = "The following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser
